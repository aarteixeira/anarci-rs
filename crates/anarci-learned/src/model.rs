//! Small non-autoregressive transformer encoder (candle) that predicts the per-residue
//! label distribution. Feature-gated behind `model` so the default build carries no ML
//! runtime. The forward pass is a single parallel encode (no decoder) — the source of the
//! throughput advantage over ANARCII's autoregressive decode.
//!
//! This is the pilot architecture: char embedding + learned positions + N pre-norm
//! self-attention/FFN blocks + a linear head over [`crate::decode::NUM_LABELS`]. It is
//! deliberately tiny; the P3 gate is only that it *learns* (overfits a small set) and that
//! its logits, decoded, are always valid.

use crate::decode::NUM_LABELS;
use candle_core::{DType, Device, Result, Tensor, D};
use candle_nn::{
    embedding, layer_norm, linear, loss, ops, AdamW, Embedding, LayerNorm, Linear, Module,
    Optimizer, ParamsAdamW, VarBuilder, VarMap,
};

/// 0 = pad, 1..=20 = the 20 standard amino acids, 21 = other/unknown.
pub const VOCAB: usize = 22;
const AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
const MAX_LEN: usize = 1024;

/// Map a residue sequence to token ids.
pub fn tokenize(seq: &[u8]) -> Vec<u32> {
    seq.iter()
        .map(|&b| match AA.iter().position(|&a| a == b.to_ascii_uppercase()) {
            Some(i) => (i + 1) as u32,
            None => 21,
        })
        .collect()
}

struct Block {
    ln1: LayerNorm,
    qkv: Linear,
    proj: Linear,
    ln2: LayerNorm,
    ff1: Linear,
    ff2: Linear,
    dim: usize,
}

impl Block {
    fn new(dim: usize, ff: usize, vb: VarBuilder) -> Result<Self> {
        Ok(Self {
            ln1: layer_norm(dim, 1e-5, vb.pp("ln1"))?,
            qkv: linear(dim, 3 * dim, vb.pp("qkv"))?,
            proj: linear(dim, dim, vb.pp("proj"))?,
            ln2: layer_norm(dim, 1e-5, vb.pp("ln2"))?,
            ff1: linear(dim, ff, vb.pp("ff1"))?,
            ff2: linear(ff, dim, vb.pp("ff2"))?,
            dim,
        })
    }

    fn forward(&self, x: &Tensor) -> Result<Tensor> {
        // pre-norm single-head self-attention (bidirectional: no causal mask)
        let h = self.ln1.forward(x)?;
        let qkv = self.qkv.forward(&h)?; // [B,T,3D]
        // `.contiguous()` is required by the CUDA matmul backend (cuBLAS): `narrow` and
        // `transpose` yield non-contiguous views that the CPU backend tolerates but CUDA
        // rejects ("matmul is only supported for contiguous tensors").
        let q = qkv.narrow(D::Minus1, 0, self.dim)?.contiguous()?;
        let k = qkv.narrow(D::Minus1, self.dim, self.dim)?.contiguous()?;
        let v = qkv.narrow(D::Minus1, 2 * self.dim, self.dim)?.contiguous()?;
        let scale = (self.dim as f64).sqrt();
        let att = (q.matmul(&k.transpose(D::Minus2, D::Minus1)?.contiguous()?)? / scale)?; // [B,T,T]
        let att = ops::softmax(&att, D::Minus1)?;
        let ctx = att.matmul(&v)?; // [B,T,D]
        let x = (x + self.proj.forward(&ctx)?)?;
        // pre-norm FFN
        let h2 = self.ln2.forward(&x)?;
        let ff = self.ff2.forward(&self.ff1.forward(&h2)?.gelu()?)?;
        x + ff
    }
}

/// Fixed sinusoidal positional encoding `[max_len, dim]` — a strong, deterministic position
/// signal (unlike a zero-init learned table that the token embedding drowns out).
fn sinusoidal(max_len: usize, dim: usize, device: &Device) -> Result<Tensor> {
    let mut data = vec![0f32; max_len * dim];
    for pos in 0..max_len {
        for i in 0..dim {
            let div = 10000f32.powf((2 * (i / 2)) as f32 / dim as f32);
            let a = pos as f32 / div;
            data[pos * dim + i] = if i % 2 == 0 { a.sin() } else { a.cos() };
        }
    }
    Tensor::from_vec(data, (max_len, dim), device)
}

/// The encoder-only numbering model.
pub struct Model {
    emb: Embedding,
    pos: Tensor, // fixed sinusoidal positional encoding [MAX_LEN, dim]
    blocks: Vec<Block>,
    lnf: LayerNorm,
    head: Linear,
}

impl Model {
    pub fn new(dim: usize, layers: usize, ff: usize, vb: VarBuilder) -> Result<Self> {
        let blocks = (0..layers)
            .map(|i| Block::new(dim, ff, vb.pp(format!("block{i}"))))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            emb: embedding(VOCAB, dim, vb.pp("emb"))?,
            pos: sinusoidal(MAX_LEN, dim, vb.device())?,
            blocks,
            lnf: layer_norm(dim, 1e-5, vb.pp("lnf"))?,
            head: linear(dim, NUM_LABELS, vb.pp("head"))?,
        })
    }

    /// tokens `[B,T]` (u32) → label logits `[B,T,NUM_LABELS]`.
    pub fn forward(&self, tokens: &Tensor) -> Result<Tensor> {
        let (_b, t) = tokens.dims2()?;
        let mut x = self.emb.forward(tokens)?; // [B,T,D]
        let pos = self.pos.narrow(0, 0, t)?.unsqueeze(0)?; // [1,T,D]
        x = x.broadcast_add(&pos)?;
        for blk in &self.blocks {
            x = blk.forward(&x)?;
        }
        self.head.forward(&self.lnf.forward(&x)?)
    }
}

/// Per-residue logits for one sequence as a `Vec<Vec<f32>>` (for [`crate::viterbi_decode`]).
pub fn logits_rows(model: &Model, seq: &[u8], device: &Device) -> Result<Vec<Vec<f32>>> {
    let toks = Tensor::new(tokenize(seq), device)?.unsqueeze(0)?; // [1,T]
    let logits = model.forward(&toks)?.squeeze(0)?; // [T,NUM_LABELS]
    logits.to_vec2::<f32>()
}

/// Overfit a tiny labeled set (each item = (tokens, target label indices)) and return the
/// final per-residue training accuracy. The P3 learn-ability gate: this must reach ~1.0.
pub fn overfit(
    examples: &[(Vec<u32>, Vec<u32>)],
    dim: usize,
    layers: usize,
    ff: usize,
    steps: usize,
    lr: f64,
    device: &Device,
) -> Result<(Model, f32)> {
    let varmap = VarMap::new();
    let vb = VarBuilder::from_varmap(&varmap, DType::F32, device);
    let model = Model::new(dim, layers, ff, vb)?;
    let mut opt = AdamW::new(varmap.all_vars(), ParamsAdamW { lr, ..Default::default() })?;

    for step in 0..steps {
        let mut epoch_loss = 0f32;
        for (toks, tgts) in examples {
            let t = toks.len();
            let x = Tensor::new(toks.as_slice(), device)?.unsqueeze(0)?; // [1,T]
            let logits = model.forward(&x)?.reshape((t, NUM_LABELS))?; // [T,C]
            let y = Tensor::new(tgts.as_slice(), device)?; // [T]
            let l = loss::cross_entropy(&logits, &y)?;
            epoch_loss += l.to_scalar::<f32>()?;
            opt.backward_step(&l)?;
        }
        if std::env::var("OVERFIT_LOG").is_ok() && (step % 50 == 0 || step + 1 == steps) {
            eprintln!("step {step:4}  loss {:.4}", epoch_loss / examples.len() as f32);
        }
    }

    // final accuracy
    let (mut correct, mut total) = (0usize, 0usize);
    for (toks, tgts) in examples {
        let t = toks.len();
        let x = Tensor::new(toks.as_slice(), device)?.unsqueeze(0)?;
        let pred = model.forward(&x)?.reshape((t, NUM_LABELS))?.argmax(D::Minus1)?;
        let pred = pred.to_vec1::<u32>()?;
        for (p, g) in pred.iter().zip(tgts) {
            total += 1;
            if p == g {
                correct += 1;
            }
        }
    }
    let acc = correct as f32 / total.max(1) as f32;
    Ok((model, acc))
}
