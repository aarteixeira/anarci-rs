//! Full training driver for the learned numbering model (P5/P6). Trains the encoder on the
//! leakage-safe dataset, evaluates on the val split each epoch, and checkpoints the
//! best-on-val model to safetensors. Device auto-selects CUDA -> Metal -> CPU.
//!
//!   CPU smoke:  cargo run -p anarci-learned --example train --features model --release
//!   GPU (gnode2): cargo run -p anarci-learned --example train --features model,cuda --release
//!
//! Config via env: ANARCI_EPOCHS, ANARCI_DIM, ANARCI_LAYERS, ANARCI_FF, ANARCI_LR,
//! ANARCI_OUT_DIR, ANARCI_SMOKE (=1 uses a tiny subset for a fast end-to-end check).

use anarci_core::{number_sequence_from_alignment, State, StateType};
use anarci_learned::data::{build_dataset, DomainRec, Split};
use anarci_learned::decode::label_to_index;
use anarci_learned::model::{tokenize, Model};
use anarci_learned::{decode_sequence, viterbi_decode};
use candle_core::{DType, Device, Tensor, D};
use candle_nn::{loss, AdamW, Optimizer, ParamsAdamW, VarBuilder, VarMap};
use flate2::read::GzDecoder;
use serde::Deserialize;
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Deserialize)]
struct Golden {
    sequences: Vec<SeqRec>,
}
#[derive(Deserialize)]
struct SeqRec {
    seq: String,
    domains: Vec<DomRec>,
}
#[derive(Deserialize)]
struct DomRec {
    state_vector: Vec<((u8, String), Option<usize>)>,
    chain_type: String,
}

fn env<T: std::str::FromStr>(k: &str, default: T) -> T {
    std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(default)
}

fn device() -> Device {
    #[cfg(feature = "cuda")]
    if let Ok(d) = Device::new_cuda(0) {
        eprintln!("device: cuda");
        return d;
    }
    if let Ok(d) = Device::new_metal(0) {
        eprintln!("device: metal");
        return d;
    }
    eprintln!("device: cpu");
    Device::Cpu
}

/// One training/eval item: tokens, per-residue target label indices, residues, chain type.
struct Item {
    toks: Vec<u32>,
    tgts: Vec<u32>,
    residues: Vec<u8>,
    chain_type: String,
    kind: &'static str,
}

fn main() -> candle_core::Result<()> {
    let dim: usize = env("ANARCI_DIM", 256);
    let layers: usize = env("ANARCI_LAYERS", 4);
    let ff: usize = env("ANARCI_FF", 512);
    let lr: f64 = env("ANARCI_LR", 3e-4);
    let epochs: usize = env("ANARCI_EPOCHS", 40);
    let out_dir: String = env("ANARCI_OUT_DIR", "ckpt".to_string());
    let smoke: usize = env("ANARCI_SMOKE", 0);

    // --- load oracle domains + build the leakage-safe dataset --------------------------
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/golden.json.gz");
    let g: Golden = serde_json::from_reader(GzDecoder::new(std::fs::File::open(&path).unwrap())).unwrap();
    // chain type per (seq_idx, dom_idx) so we can attach it to examples after augmentation
    let mut doms = Vec::new();
    for s in &g.sequences {
        for d in &s.domains {
            let sv: Vec<State> = d
                .state_vector
                .iter()
                .map(|((id, t), si)| State {
                    id: *id,
                    typ: match t.as_str() {
                        "m" => StateType::M,
                        "i" => StateType::I,
                        _ => StateType::D,
                    },
                    si: *si,
                })
                .collect();
            doms.push((DomainRec { seq: s.seq.as_bytes().to_vec(), sv }, d.chain_type.clone()));
        }
    }
    let domain_recs: Vec<DomainRec> = doms.iter().map(|(d, _)| d.clone()).collect();
    let examples = build_dataset(&domain_recs, &[(66, 128), (39, 128)]);

    // convert to Items; use a representative chain type (H) for multi/junk (only used for
    // scheme-specific numbering during the validity check, not for the label targets).
    let to_items = |split: Split| -> Vec<Item> {
        examples
            .iter()
            .filter(|e| e.split == split)
            .map(|e| Item {
                toks: tokenize(&e.data.seq),
                tgts: e.data.labels.iter().map(|&l| label_to_index(l) as u32).collect(),
                residues: e.data.seq.clone(),
                chain_type: "H".to_string(),
                kind: e.kind,
            })
            .collect()
    };
    let mut train = to_items(Split::Train);
    let val = to_items(Split::Val);
    if smoke > 0 {
        train.truncate(16);
    }
    eprintln!("train {} / val {} examples  (smoke={})", train.len(), val.len(), smoke);
    // Fixed sample of clean-domain TRAIN examples: measuring single_acc here vs on val is the
    // underfit-vs-overfit diagnostic behind the go/no-go (both low ⇒ underfit/architecture;
    // train ≫ val ⇒ data-limited, which P6's data scaling addresses).
    let train_probe: Vec<&Item> = train.iter().filter(|it| it.kind == "single").take(200).collect();

    // --- model + optimizer --------------------------------------------------------------
    let dev = device();
    let varmap = VarMap::new();
    let vb = VarBuilder::from_varmap(&varmap, DType::F32, &dev);
    let model = Model::new(dim, layers, ff, vb)?;
    let mut opt = AdamW::new(varmap.all_vars(), ParamsAdamW { lr, ..Default::default() })?;
    std::fs::create_dir_all(&out_dir).ok();

    let n_labels = anarci_learned::decode::NUM_LABELS;
    let mut best_val = -1.0f32;
    let epochs = if smoke > 0 { epochs.min(3) } else { epochs };
    // Minibatch via gradient accumulation: stack `accum` per-example losses into one
    // backward step (candle backprops through all of them), giving a smoothed gradient over
    // variable-length sequences with no padding/masking. accum=1 == the original B=1 loop.
    let accum: usize = std::env::var("ANARCI_ACCUM").ok().and_then(|s| s.parse().ok()).unwrap_or(1).max(1);
    eprintln!("grad-accum (minibatch) = {accum}");

    for ep in 0..epochs {
        let mut loss_sum = 0f32;
        let mut mb: Vec<Tensor> = Vec::with_capacity(accum);
        for it in &train {
            let t = it.toks.len();
            let x = Tensor::new(it.toks.as_slice(), &dev)?.unsqueeze(0)?;
            let logits = model.forward(&x)?.reshape((t, n_labels))?;
            let y = Tensor::new(it.tgts.as_slice(), &dev)?;
            let l = loss::cross_entropy(&logits, &y)?;
            loss_sum += l.to_scalar::<f32>()?;
            mb.push(l);
            if mb.len() == accum {
                opt.backward_step(&Tensor::stack(&mb, 0)?.mean_all()?)?;
                mb.clear();
            }
        }
        if !mb.is_empty() {
            opt.backward_step(&Tensor::stack(&mb, 0)?.mean_all()?)?;
            mb.clear();
        }

        // --- val: per-residue accuracy, EXACT numbering match, and decoded validity,
        // broken down by example kind. The go/no-go metric is the `single` kind (clean
        // domains): `exact` there is the fraction numbered byte-identically to the ANARCI
        // oracle — i.e. parity. Slots per kind: [res_correct, res_total, exact, valid, n].
        let mut by: HashMap<&str, [usize; 5]> = HashMap::new();
        for it in &val {
            let t = it.toks.len();
            let x = Tensor::new(it.toks.as_slice(), &dev)?.unsqueeze(0)?;
            let logits = model.forward(&x)?.reshape((t, n_labels))?;
            let pred = logits.argmax(D::Minus1)?.to_vec1::<u32>()?;
            let e = by.entry(it.kind).or_insert([0; 5]);
            for (p, gt) in pred.iter().zip(&it.tgts) {
                e[1] += 1;
                if p == gt {
                    e[0] += 1;
                }
            }
            // decoded labels = what the engine actually emits: exact match + validity
            let rows = logits.to_vec2::<f32>()?;
            let labels = viterbi_decode(&rows);
            if labels.iter().zip(&it.tgts).all(|(l, gt)| label_to_index(*l) as u32 == *gt) {
                e[2] += 1;
            }
            let ok = decode_sequence(&labels).iter().all(|sv| {
                let ms: Vec<u8> = sv.iter().filter(|s| s.typ == StateType::M).map(|s| s.id).collect();
                ms.windows(2).all(|w| w[0] < w[1])
                    && number_sequence_from_alignment(sv, &it.residues, "imgt", Some(&it.chain_type)).is_ok()
            });
            if ok {
                e[3] += 1;
            }
            e[4] += 1;
        }
        let s = |slot: usize| by.values().map(|v| v[slot]).sum::<usize>();
        let acc = s(0) as f32 / s(1).max(1) as f32;
        let valrate = s(3) as f32 / s(4).max(1) as f32;
        let sg = by.get("single").copied().unwrap_or([0; 5]);
        let single_acc = sg[0] as f32 / sg[1].max(1) as f32;
        let single_exact = sg[2] as f32 / sg[4].max(1) as f32;
        // train single_acc on the fixed probe (same clean-domain metric) — the gap to val
        // single_acc is the underfit/overfit diagnostic.
        let (mut tc, mut tt) = (0usize, 0usize);
        for it in &train_probe {
            let t = it.toks.len();
            let x = Tensor::new(it.toks.as_slice(), &dev)?.unsqueeze(0)?;
            let pred = model.forward(&x)?.reshape((t, n_labels))?.argmax(D::Minus1)?.to_vec1::<u32>()?;
            for (p, gt) in pred.iter().zip(&it.tgts) {
                tt += 1;
                if *p == *gt {
                    tc += 1;
                }
            }
        }
        let train_single_acc = tc as f32 / tt.max(1) as f32;
        eprintln!(
            "epoch {ep:3}  train_loss {:.4}  val_acc {:.4}  val_valid {:.4}  single_acc {:.4}  single_exact {:.4}  train_single_acc {:.4}",
            loss_sum / train.len().max(1) as f32,
            acc,
            valrate,
            single_acc,
            single_exact,
            train_single_acc
        );
        if ep % 10 == 0 || ep + 1 == epochs {
            let mut kinds: Vec<&str> = by.keys().copied().collect();
            kinds.sort();
            for k in kinds {
                let v = by[k];
                eprintln!(
                    "    [{k:8}] n={:4}  res_acc {:.4}  exact {:.4}  valid {:.4}",
                    v[4],
                    v[0] as f32 / v[1].max(1) as f32,
                    v[2] as f32 / v[4].max(1) as f32,
                    v[3] as f32 / v[4].max(1) as f32,
                );
            }
        }
        // checkpoint on best overall per-residue accuracy (unchanged criterion)
        if acc > best_val {
            best_val = acc;
            varmap.save(format!("{out_dir}/best.safetensors"))?;
        }
    }
    eprintln!("done. best val_acc {best_val:.4}  checkpoint {out_dir}/best.safetensors");
    Ok(())
}
