//! Shared sequence, structure, annotation and CLI entry points.
pub mod cif;
pub mod engine;
pub mod output;
pub mod pairing;
pub mod pdb;
pub mod seqio;
pub mod structure;
pub mod translate;

use anarci_core::orchestrate::{default_allow, run_anarci, GermlineMethod, HmmEngine, SeqResult};
use clap::{Parser, ValueEnum};
use output::DomainOut;
use rayon::prelude::*;
use seqio::Record;
use std::collections::BTreeSet;
use std::io::{Read, Write};
use translate::FrameOrder;

#[derive(Debug, Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum Eng {
    /// FEW.hmm — pan-species, fast (default; right for NGS).
    Pan,
    /// ALL.hmm — every species, byte-exact ANARCI parity.
    All,
}
#[derive(Debug, Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum Fmt {
    Tsv,
    Json,
    /// ANARCI-style per-position table.
    Anarci,
}
#[derive(Debug, Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum InType {
    Auto,
    Protein,
    Dna,
    Pdb,
    /// mmCIF coordinates (also accepts mmcif).
    #[value(alias = "mmcif")]
    Cif,
}

/// Annotate antibody/TCR V domains from protein/DNA/PDB inputs.
#[derive(Debug, Clone, Parser)]
#[command(name = "anarci", version, about)]
pub struct Args {
    /// Input file (FASTA/FASTQ/PDB). Omit or "-" for stdin.
    pub input: Option<String>,
    /// Numbering scheme: imgt, kabat, chothia, martin, aho.
    #[arg(long, default_value = "imgt")]
    pub scheme: String,
    /// HMM engine: pan (fast) or all (exact).
    #[arg(long, value_enum, default_value_t = Eng::Pan)]
    pub engine: Eng,
    /// Output format.
    #[arg(long, value_enum, default_value_t = Fmt::Tsv)]
    pub format: Fmt,
    /// Output file (default: stdout).
    #[arg(short, long)]
    pub output: Option<String>,
    /// Force input interpretation (default: auto-detect by extension + content).
    #[arg(long, value_enum, default_value_t = InType::Auto)]
    pub input_type: InType,
    /// HMM bit-score threshold. Lower it (e.g. 20) to recover truncated/partial V.
    #[arg(long, default_value_t = 80.0)]
    pub bit_score_threshold: f64,
    /// Assign V/J germline genes.
    #[arg(long)]
    pub assign_germline: bool,
    /// Restrict to these chain types (comma-separated: H,K,L,A,B,G,D). Default: all.
    #[arg(long, value_delimiter = ',')]
    pub restrict: Vec<String>,
    /// Direct API chain filter; an empty set selects no chains.
    #[arg(skip)]
    pub allow: Option<BTreeSet<String>>,
    /// After DNA input, print reading-frame hit counts to stderr.
    #[arg(long)]
    pub report_frames: bool,
    /// For PDB/mmCIF input: infer VH/VL pairing and antibody format (scFv, diabody,
    /// tandem, VHH, mixtures) from 3D Cα contacts, instead of numbering.
    #[arg(long)]
    pub pair: bool,
    /// Restrict species assignment (comma-separated), reverting to any species when none
    /// scores well (ANARCI's behavior — avoids over-calling camelid on conventional chains).
    /// "all" disables limiting.
    #[arg(long, default_value = "human,mouse")]
    pub species: String,
    /// Germline method: identity or evalue (default: evalue for pan, identity for all).
    #[arg(long)]
    pub germline_method: Option<String>,
    /// Worker threads; 0 uses all available Rayon workers.
    #[arg(long, default_value_t = 0)]
    pub ncpu: usize,
}

impl Default for Args {
    fn default() -> Self {
        Self::parse_from(["anarci"])
    }
}

impl Args {
    pub fn validate(&self) -> Result<(), String> {
        anarci_core::orchestrate::resolve_scheme(&self.scheme).map_err(|e| e.to_string())?;
        self.method()?;
        if !self.bit_score_threshold.is_finite() || self.bit_score_threshold < 0.0 {
            return Err("bit_score_threshold must be finite and non-negative".into());
        }
        for c in self.restrict.iter().chain(self.allow.iter().flatten()) {
            if !default_allow().contains(&c.to_ascii_uppercase()) {
                return Err(format!("unknown chain type: {c}"));
            }
        }
        if self.species.split(',').any(|s| s.trim().is_empty()) {
            return Err("species must contain non-empty names or 'all'".into());
        }
        Ok(())
    }
    pub fn method(&self) -> Result<GermlineMethod, String> {
        match self.germline_method.as_deref() {
            Some(m) => GermlineMethod::parse(m).map_err(|e| e.to_string()),
            None => Ok(if self.engine == Eng::Pan {
                GermlineMethod::Evalue
            } else {
                GermlineMethod::Identity
            }),
        }
    }
    pub fn allowed_species(&self) -> Option<Vec<String>> {
        if self.species.eq_ignore_ascii_case("all") {
            None
        } else {
            Some(
                self.species
                    .split(',')
                    .map(|s| s.trim().to_lowercase())
                    .collect(),
            )
        }
    }
    fn allow(&self) -> BTreeSet<String> {
        if let Some(allow) = &self.allow {
            return allow.iter().map(|s| s.to_ascii_uppercase()).collect();
        }
        if self.restrict.is_empty() {
            default_allow()
        } else {
            self.restrict
                .iter()
                .map(|s| s.to_ascii_uppercase())
                .collect()
        }
    }
}

fn with_workers<T: Send>(
    args: &Args,
    run: impl FnOnce() -> Result<T, String> + Send,
) -> Result<T, String> {
    if args.ncpu == 0 {
        run()
    } else {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(args.ncpu)
            .build()
            .map_err(|e| format!("create worker pool: {e}"))?;
        pool.install(run)
    }
}

pub fn parse_records(raw: &[u8]) -> Result<Vec<Record>, String> {
    seqio::read_records(std::io::BufReader::new(raw))
}

fn validate_records(records: &[Record], dna: bool) -> Result<(), String> {
    let alphabet: &[u8] = if dna {
        b"ACGTUNRYKMSWBDHV"
    } else {
        b"ACDEFGHIKLMNPQRSTVWYBXZJUO*"
    };
    for r in records {
        if r.id.trim().is_empty() {
            return Err("record identifier must not be empty".into());
        }
        if let Some(b) = r
            .seq
            .iter()
            .find(|b| !alphabet.contains(&b.to_ascii_uppercase()))
        {
            return Err(format!(
                "record {}: invalid {} character {:?}",
                r.id,
                if dna { "DNA" } else { "protein" },
                *b as char
            ));
        }
    }
    Ok(())
}

fn engine_for(args: &Args) -> Result<&'static anarci_hmm::Engine, String> {
    engine::load(if args.engine == Eng::Pan {
        engine::EngineKind::Pan
    } else {
        engine::EngineKind::All
    })
}

fn annotate_with(
    eng: &(dyn HmmEngine + Sync),
    records: &[Record],
    args: &Args,
    scheme: &str,
    germline: bool,
) -> Result<Vec<SeqResult>, String> {
    let batch: Vec<_> = records
        .iter()
        .map(|r| (r.name(), r.seq.to_ascii_uppercase()))
        .collect();
    run_anarci(
        eng,
        &batch,
        scheme,
        &args.allow(),
        germline,
        args.allowed_species().as_deref(),
        args.bit_score_threshold,
        args.engine == Eng::Pan,
        args.method()?,
    )
    .map_err(|e| format!("numbering failed: {e}"))
}

pub fn annotate_protein(records: &[Record], args: &Args) -> Result<Vec<DomainOut>, String> {
    args.validate()?;
    validate_records(records, false)?;
    if records.is_empty() {
        return Ok(Vec::new());
    }
    let eng = engine_for(args)?;
    with_workers(args, || {
        let results = annotate_with(
            &engine::Adapter(eng),
            records,
            args,
            &args.scheme,
            args.assign_germline,
        )?;
        Ok(records
            .iter()
            .zip(&results)
            .flat_map(|(r, res)| output::rows(&r.name(), res, None))
            .collect())
    })
}

#[derive(serde::Serialize)]
pub struct DnaAnnotation {
    pub domains: Vec<DomainOut>,
    pub frame_counts: [u64; 6],
}

pub fn annotate_dna(records: &[Record], args: &Args) -> Result<DnaAnnotation, String> {
    args.validate()?;
    validate_records(records, true)?;
    let frames = FrameOrder::default();
    if records.is_empty() {
        return Ok(DnaAnnotation {
            domains: Vec::new(),
            frame_counts: frames.counts(),
        });
    }
    let eng = engine_for(args)?;
    let domains = with_workers(args, || {
        let batches: Result<Vec<Vec<DomainOut>>, String> = records
            .par_iter()
            .map(|r| number_read(&engine::Adapter(eng), r, args, &frames))
            .collect();
        Ok(batches?.into_iter().flatten().collect())
    })?;
    Ok(DnaAnnotation {
        domains,
        frame_counts: frames.counts(),
    })
}

fn number_read(
    eng: &(dyn HmmEngine + Sync),
    rec: &Record,
    args: &Args,
    frames: &FrameOrder,
) -> Result<Vec<DomainOut>, String> {
    for f in frames.order() {
        let prot: Vec<u8> = translate::translate_frame(&rec.seq, f)
            .into_iter()
            .map(|b| if b == b'*' { b'X' } else { b })
            .collect();
        if prot.len() < 10 {
            continue;
        }
        let record = Record {
            id: rec.id.clone(),
            desc: rec.desc.clone(),
            seq: prot,
        };
        let res = annotate_with(eng, &[record], args, &args.scheme, args.assign_germline)?;
        let rows = output::rows(&rec.name(), &res[0], Some(f));
        if !rows.is_empty() {
            frames.record(f);
            return Ok(rows);
        }
    }
    Ok(Vec::new())
}

/// Parse C-alpha residues from the first model. Supported formats are pdb, cif, mmcif.
pub fn parse_structure(raw: &[u8], format: &str) -> Result<Vec<structure::StructChain>, String> {
    match format.to_ascii_lowercase().as_str() {
        "pdb" => pdb::chains_coords(raw),
        "cif" | "mmcif" => cif::chains_coords(raw),
        _ => Err(format!("unknown structure format: {format}")),
    }
}

pub fn pair_structure(raw: &[u8], format: &str, args: &Args) -> Result<pairing::Pairing, String> {
    let chains = parse_structure(raw, format)?;
    pair_chains(&chains, args)
}

pub fn pair_chains(
    chains: &[structure::StructChain],
    args: &Args,
) -> Result<pairing::Pairing, String> {
    args.validate()?;
    let records: Vec<_> = chains.iter().map(|c| c.to_record()).collect();
    if records.is_empty() {
        return Ok(pairing::infer(&[]));
    }
    let eng = engine_for(args)?;
    with_workers(args, || {
        // Contact masks always use IMGT, independent of the requested display scheme.
        let results = annotate_with(&engine::Adapter(eng), &records, args, "imgt", true)?;
        let vdoms = domains_with_coords(chains, &results);
        let mut result = pairing::infer(&vdoms);
        pairing::assign_antigens(&mut result, &vdoms, chains);
        Ok(result)
    })
}

fn domains_with_coords(
    chains: &[structure::StructChain],
    results: &[SeqResult],
) -> Vec<pairing::VDomain> {
    let mut vdoms = Vec::new();
    for (chain, res) in chains.iter().zip(results) {
        let (Some(details), Some(numbered)) = (&res.details, &res.numbered) else {
            continue;
        };
        for (di, (d, nb)) in details.iter().zip(numbered).enumerate() {
            let end = d.query_end.min(chain.residues.len());
            let start = d.query_start.unwrap_or(0).min(end);
            let mut ri = nb.start.unwrap_or(start);
            let mut paratope = Vec::new();
            for r in &nb.residues {
                if r.aa == b'-' {
                    continue;
                }
                if ri >= end {
                    break;
                }
                if (27..=38).contains(&r.num)
                    || (56..=65).contains(&r.num)
                    || (105..=117).contains(&r.num)
                {
                    paratope.push(chain.residues[ri].ca);
                }
                ri += 1;
            }
            vdoms.push(pairing::VDomain {
                id: format!("{}:{}", chain.id, di),
                chain: chain.id.clone(),
                chain_type: d.chain_type.clone(),
                species: d.species.clone(),
                v_gene: d
                    .germlines
                    .as_ref()
                    .and_then(|g| g.v_gene.as_ref())
                    .map(|(_, g)| g.clone()),
                j_gene: d
                    .germlines
                    .as_ref()
                    .and_then(|g| g.j_gene.as_ref())
                    .map(|(_, g)| g.clone()),
                cas: chain.residues[start..end].iter().map(|r| r.ca).collect(),
                paratope,
            });
        }
    }
    vdoms
}

pub fn render_output(rows: &[DomainOut], format: Fmt) -> Result<String, String> {
    match format {
        Fmt::Json => serde_json::to_string_pretty(rows)
            .map(|s| s + "\n")
            .map_err(|e| e.to_string()),
        Fmt::Tsv => Ok(std::iter::once(output::tsv_header().to_string())
            .chain(rows.iter().map(output::tsv_row))
            .collect::<Vec<_>>()
            .join("\n")
            + "\n"),
        Fmt::Anarci => Ok(rows.iter().map(output::anarci_block).collect()),
    }
}

pub fn render_pairing(result: &pairing::Pairing, format: Fmt) -> Result<String, String> {
    if format == Fmt::Json {
        serde_json::to_string_pretty(result)
            .map(|s| s + "\n")
            .map_err(|e| e.to_string())
    } else {
        Ok(pairing::report(result))
    }
}

fn read_input(path: Option<&str>) -> Result<Vec<u8>, String> {
    let mut buf = Vec::new();
    match path {
        None | Some("-") => std::io::stdin()
            .read_to_end(&mut buf)
            .map_err(|e| e.to_string())?,
        Some(p) => std::fs::File::open(p)
            .map_err(|e| format!("open {p}: {e}"))?
            .read_to_end(&mut buf)
            .map_err(|e| e.to_string())?,
    };
    Ok(buf)
}

fn render_input(raw: &[u8], args: &Args) -> Result<String, String> {
    args.validate()?;
    let lower = args.input.as_deref().unwrap_or("").to_ascii_lowercase();
    let structure_format = match args.input_type {
        InType::Pdb => Some("pdb"),
        InType::Cif => Some("cif"),
        InType::Auto if lower.ends_with(".cif") || lower.ends_with(".mmcif") => Some("cif"),
        InType::Auto if lower.ends_with(".pdb") || lower.ends_with(".ent") => Some("pdb"),
        _ => None,
    };
    if args.pair && structure_format.is_none() {
        return Err("--pair requires a PDB/mmCIF structure input".into());
    }
    if let Some(format) = structure_format {
        let chains = parse_structure(raw, format)?;
        if args.pair {
            render_pairing(&pair_chains(&chains, args)?, args.format)
        } else {
            render_output(
                &annotate_protein(
                    &chains.iter().map(|c| c.to_record()).collect::<Vec<_>>(),
                    args,
                )?,
                args.format,
            )
        }
    } else {
        let records = parse_records(raw)?;
        let dna = args.input_type == InType::Dna
            || (args.input_type == InType::Auto
                && !records.is_empty()
                && records.iter().all(|r| seqio::is_nucleotide(&r.seq)));
        if dna {
            let result = annotate_dna(&records, args)?;
            if args.report_frames {
                eprintln!(
                    "frame hit counts (0-2 fwd, 3-5 rev): {:?}",
                    result.frame_counts
                );
            }
            render_output(&result.domains, args.format)
        } else {
            render_output(&annotate_protein(&records, args)?, args.format)
        }
    }
}

pub fn run(args: Args) -> Result<(), String> {
    args.validate()?;
    let raw = read_input(args.input.as_deref())?;
    let rendered = render_input(&raw, &args)?;
    match args.output {
        Some(p) => std::fs::write(&p, rendered).map_err(|e| format!("write {p}: {e}")),
        None => std::io::stdout()
            .lock()
            .write_all(rendered.as_bytes())
            .map_err(|e| e.to_string()),
    }
}

/// CLI arguments exclude the executable name. Help/version return successfully.
pub fn cli(arguments: Vec<String>) -> Result<(), String> {
    let parsed = Args::try_parse_from(std::iter::once("anarci".to_string()).chain(arguments));
    match parsed {
        Ok(args) => run(args),
        Err(e)
            if matches!(
                e.kind(),
                clap::error::ErrorKind::DisplayHelp | clap::error::ErrorKind::DisplayVersion
            ) =>
        {
            e.print().map_err(|e| e.to_string())
        }
        Err(e) => Err(e.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anarci_core::{CResult, CoreError, Hsp};

    struct Broken;
    impl HmmEngine for Broken {
        fn scan_one(&self, _: &str, _: &[u8]) -> Vec<Hsp> {
            panic!("checked scan required")
        }
        fn scan_one_checked(&self, _: &str, _: &[u8]) -> CResult<Vec<Hsp>> {
            Err(CoreError::Assertion("synthetic scan error".into()))
        }
    }
    struct Synthetic {
        chain: &'static str,
    }
    impl HmmEngine for Synthetic {
        fn scan_one(&self, _: &str, seq: &[u8]) -> Vec<Hsp> {
            vec![Hsp {
                hit_id: format!("human_{}", self.chain),
                hit_description: String::new(),
                evalue: 1e-20,
                bitscore: 200.,
                bias: 0.,
                query_start: 0,
                query_end: seq.len(),
                hit_start: 0,
                hit_end: seq.len(),
                rf: "x".repeat(seq.len()),
                pp: "*".repeat(seq.len()),
                order: 0,
            }]
        }
    }
    fn synthetic_record() -> Record {
        Record {
            id: "synthetic".into(),
            desc: String::new(),
            seq: vec![b'A'; 128],
        }
    }

    #[test]
    fn invalid_config_and_inputs_fail_before_compute() {
        for scheme in ["unknown", ""] {
            let args = Args {
                scheme: scheme.into(),
                ..Args::default()
            };
            assert!(annotate_protein(&[], &args).is_err());
            assert!(annotate_dna(&[], &args).is_err());
        }
        for threshold in [f64::NAN, f64::INFINITY, -1.] {
            assert!(Args {
                bit_score_threshold: threshold,
                ..Args::default()
            }
            .validate()
            .is_err());
        }
        let args = Args::default();
        let invalid = Record {
            id: "invalid".into(),
            desc: String::new(),
            seq: b"AC?".to_vec(),
        };
        assert!(annotate_protein(&[invalid.clone()], &args).is_err());
        assert!(annotate_dna(&[invalid], &args).is_err());
        assert_eq!(args.method().unwrap(), GermlineMethod::Evalue);
        assert_eq!(
            Args {
                engine: Eng::All,
                ..args.clone()
            }
            .method()
            .unwrap(),
            GermlineMethod::Identity
        );
        assert!(Args {
            germline_method: Some("bad".into()),
            ..args
        }
        .validate()
        .is_err());
    }
    #[test]
    fn propagates_protein_and_dna_scan_errors() {
        let args = Args {
            engine: Eng::All,
            ..Args::default()
        };
        assert!(
            annotate_with(&Broken, &[synthetic_record()], &args, "imgt", false)
                .unwrap_err()
                .contains("synthetic scan error")
        );
        let record = Record {
            id: "dna".into(),
            desc: String::new(),
            seq: b"GCT".repeat(128),
        };
        assert!(number_read(&Broken, &record, &args, &FrameOrder::default())
            .err()
            .unwrap()
            .contains("synthetic scan error"));
    }
    #[test]
    fn synthetic_numbering_all_schemes_and_chain_classes() {
        let args = Args {
            engine: Eng::All,
            species: "all".into(),
            ..Args::default()
        };
        for scheme in ["imgt", "kabat", "chothia", "martin", "aho"] {
            let results = annotate_with(
                &Synthetic { chain: "H" },
                &[synthetic_record()],
                &args,
                scheme,
                false,
            )
            .unwrap();
            let rows = output::rows("synthetic", &results[0], None);
            assert_eq!(rows.len(), 1);
            assert!(render_output(&rows, Fmt::Json)
                .unwrap()
                .contains("synthetic"));
            assert!(render_output(&rows, Fmt::Tsv)
                .unwrap()
                .starts_with("query\tdomain"));
            assert!(render_output(&rows, Fmt::Anarci).unwrap().ends_with("//\n"));
        }
        let no_chains = Args {
            allow: Some(BTreeSet::new()),
            ..args.clone()
        };
        let res = annotate_with(
            &Synthetic { chain: "H" },
            &[synthetic_record()],
            &no_chains,
            "imgt",
            false,
        )
        .unwrap();
        assert!(res[0].numbered.is_none());
        for chain in ["H", "K", "L", "A", "B", "G", "D"] {
            let res = annotate_with(
                &Synthetic { chain },
                &[synthetic_record()],
                &args,
                "imgt",
                false,
            )
            .unwrap();
            assert_eq!(res[0].details.as_ref().unwrap()[0].chain_type, chain);
        }
    }
    #[test]
    fn imgt_contact_mask_uses_residue_coordinates() {
        let args = Args {
            engine: Eng::All,
            ..Args::default()
        };
        let res = annotate_with(
            &Synthetic { chain: "H" },
            &[synthetic_record()],
            &args,
            "imgt",
            false,
        )
        .unwrap();
        let chains = vec![structure::StructChain {
            id: "synthetic".into(),
            residues: (0..128)
                .map(|i| structure::ResCoord {
                    aa: b'A',
                    ca: [i as f32, 0., 0.],
                })
                .collect(),
        }];
        let domains = domains_with_coords(&chains, &res);
        assert_eq!(domains[0].paratope.len(), 35);
        assert_eq!(domains[0].paratope[0], [26., 0., 0.]);
        assert_eq!(domains[0].paratope[34], [116., 0., 0.]);
    }
    #[test]
    fn pairing_output_honors_output_path() {
        let path =
            std::env::temp_dir().join(format!("anarci_pair_output_{}.json", std::process::id()));
        let input = path.with_extension("pdb");
        std::fs::write(&input, "END\n").unwrap();
        run(Args {
            input: Some(input.to_string_lossy().into()),
            input_type: InType::Pdb,
            pair: true,
            format: Fmt::Json,
            output: Some(path.to_string_lossy().into()),
            ..Args::default()
        })
        .unwrap();
        let result: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
        assert_eq!(result["n_domains"], 0);
        std::fs::remove_file(path).unwrap();
        std::fs::remove_file(input).unwrap();
    }
    #[test]
    fn actual_engines_accept_synthetic_no_hits_and_worker_counts() {
        let record = Record {
            id: "synthetic".into(),
            desc: String::new(),
            seq: b"ACDEFGHIKLMNPQ".to_vec(),
        };
        for engine in [Eng::Pan, Eng::All] {
            for ncpu in [0, 1, 2] {
                let args = Args {
                    engine,
                    ncpu,
                    ..Args::default()
                };
                assert!(annotate_protein(&[record.clone()], &args)
                    .unwrap()
                    .is_empty());
            }
        }
    }
    #[test]
    fn explicit_structure_type_overrides_suffix_and_supports_stdin_cif() {
        let args = Args {
            input: Some("synthetic.cif".into()),
            input_type: InType::Pdb,
            format: Fmt::Json,
            ..Args::default()
        };
        assert_eq!(render_input(b"END\n", &args).unwrap(), "[]\n");
        assert!(render_input(
            b"END\n",
            &Args {
                input_type: InType::Auto,
                ..args.clone()
            }
        )
        .is_err());
        assert!(render_input(
            b"END\n",
            &Args {
                input: Some("synthetic.pdb".into()),
                input_type: InType::Cif,
                ..args.clone()
            }
        )
        .is_err());
        let cif = b"data_synthetic\nloop_\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\nCA ALA A 1\n";
        assert_eq!(
            render_input(
                cif,
                &Args {
                    input: Some("-".into()),
                    input_type: InType::Cif,
                    ..args
                }
            )
            .unwrap(),
            "[]\n"
        );
        for name in ["cif", "mmcif"] {
            assert_eq!(
                Args::try_parse_from(["anarci", "--input-type", name])
                    .unwrap()
                    .input_type,
                InType::Cif
            );
        }
    }
    #[test]
    fn tsv_quotes_tabs_newlines_and_embedded_quotes() {
        let args = Args {
            engine: Eng::All,
            ..Args::default()
        };
        let results = annotate_with(
            &Synthetic { chain: "H" },
            &[synthetic_record()],
            &args,
            "imgt",
            false,
        )
        .unwrap();
        let mut rows = output::rows("id\t\"quote\"\nline", &results[0], None);
        rows[0].v_gene = Some("v\rname".into());
        rows[0].j_gene = Some("j\"name".into());
        rows[0].numbered_seq = "AA\tQQ".into();
        assert_eq!(output::tsv_row(&rows[0]), "\"id\t\"\"quote\"\"\nline\"\t0\tH\thuman\timgt\t0\t128\tfull\t-\t\"v\rname\"\t\"j\"\"name\"\t200.0\t\"AA\tQQ\"");
    }
    #[test]
    fn complete_flag_preserves_partial_boundary_status() {
        let args = Args {
            engine: Eng::All,
            ..Args::default()
        };
        let record = Record {
            seq: vec![b'A'; 127],
            ..synthetic_record()
        };
        let result =
            annotate_with(&Synthetic { chain: "H" }, &[record], &args, "imgt", false).unwrap();
        let rows = output::rows("synthetic", &result[0], None);
        assert!(rows[0].complete);
        assert_eq!(rows[0].regions.last(), Some(&("fr4", "partial")));
    }
}
