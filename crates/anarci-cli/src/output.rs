//! Output formatting: one logical row per numbered domain. JSON carries the full
//! per-position numbering; TSV is a one-line-per-domain summary; `anarci` is the
//! classic per-position table.

use anarci_core::orchestrate::{DomainInfo, SeqResult};
use anarci_core::regions::IMGT_REGIONS;
use anarci_core::Numbered;
use serde::Serialize;

/// A numbered domain, flattened for output. `frame` is set only for DNA input.
#[derive(Serialize)]
pub struct DomainOut {
    pub query: String,
    pub domain: usize,
    pub chain_type: String,
    pub species: String,
    pub scheme: String,
    /// Domain location within the input polypeptide (0-based residue indices).
    pub start: Option<usize>,
    pub end: usize,
    pub evalue: f64,
    pub bitscore: f64,
    /// True when all seven IMGT regions are represented at least partly.
    /// Inspect `regions` for partial coverage and missing terminal residues.
    pub complete: bool,
    /// Per-region completeness: "complete" | "partial" | "absent" (FR1..FR4).
    pub regions: Vec<(&'static str, &'static str)>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub frame: Option<usize>,
    // germline genes (present only with --assign-germline)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub v_gene: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub v_identity: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub j_gene: Option<String>,
    /// The domain's residues in scheme order (gaps dropped): scheme label + amino acid.
    pub numbering: Vec<Position>,
    /// The domain amino-acid sequence (gaps dropped) — handy for TSV / grep.
    pub numbered_seq: String,
    /// Every scheme position incl. gaps ('-'), for the ANARCI table. Not serialized.
    #[serde(skip)]
    pub all_positions: Vec<(String, char)>,
}

#[derive(Serialize)]
pub struct Position {
    /// Scheme number plus insertion code, e.g. "111" or "111A".
    pub n: String,
    pub aa: char,
}

/// Build the output rows for one input sequence's result.
pub fn rows(query: &str, res: &SeqResult, frame: Option<usize>) -> Vec<DomainOut> {
    let (Some(details), Some(numbered)) = (res.details.as_ref(), res.numbered.as_ref()) else {
        return Vec::new();
    };
    details
        .iter()
        .zip(numbered.iter())
        .enumerate()
        .map(|(i, (d, n))| domain_row(query, i, d, n, frame))
        .collect()
}

fn domain_row(
    query: &str,
    i: usize,
    d: &DomainInfo,
    n: &Numbered,
    frame: Option<usize>,
) -> DomainOut {
    let statuses = &d.regions.statuses;
    let regions: Vec<(&'static str, &'static str)> = IMGT_REGIONS
        .iter()
        .zip(statuses.iter())
        .map(|(&(name, _, _), s)| (name, s.as_str()))
        .collect();
    // Keep the established "full" TSV label for domains with all seven regions represented.
    // Individual region statuses retain partial boundary coverage.
    let complete = statuses.iter().all(|s| s.as_str() != "absent");

    let mut numbering = Vec::new();
    let mut numbered_seq = String::new();
    let mut all_positions = Vec::new();
    for r in &n.residues {
        let label = format!("{}{}", r.num, r.ins.trim());
        all_positions.push((label.clone(), r.aa as char));
        if r.aa != b'-' {
            numbering.push(Position {
                n: label,
                aa: r.aa as char,
            });
            numbered_seq.push(r.aa as char);
        }
    }

    let (v_gene, v_identity, j_gene) = match d.germlines.as_ref() {
        Some(g) => (
            g.v_gene.as_ref().map(|(_sp, gene)| gene.clone()),
            g.v_identity,
            g.j_gene.as_ref().map(|(_sp, gene)| gene.clone()),
        ),
        None => (None, None, None),
    };

    DomainOut {
        query: query.to_string(),
        domain: i,
        chain_type: d.chain_type.clone(),
        species: d.species.clone(),
        scheme: d.scheme.clone(),
        start: d.query_start,
        end: d.query_end,
        evalue: d.evalue,
        bitscore: d.bitscore,
        complete,
        regions,
        frame,
        v_gene,
        v_identity,
        j_gene,
        numbering,
        numbered_seq,
        all_positions,
    }
}

pub fn tsv_header() -> &'static str {
    "query\tdomain\tchain_type\tspecies\tscheme\tstart\tend\tcomplete\tframe\tv_gene\tj_gene\tbitscore\tnumbered_seq"
}

fn tsv_field(value: &str) -> String {
    if value
        .bytes()
        .any(|b| matches!(b, b'\t' | b'\r' | b'\n' | b'"'))
    {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

pub fn tsv_row(d: &DomainOut) -> String {
    let dash = || "-".to_string();
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}",
        tsv_field(&d.query),
        d.domain,
        tsv_field(&d.chain_type),
        tsv_field(&d.species),
        tsv_field(&d.scheme),
        d.start.map(|s| s.to_string()).unwrap_or_else(dash),
        d.end,
        if d.complete { "full" } else { "partial" },
        d.frame.map(|f| f.to_string()).unwrap_or_else(dash),
        d.v_gene.as_deref().map(tsv_field).unwrap_or_else(dash),
        d.j_gene.as_deref().map(tsv_field).unwrap_or_else(dash),
        d.bitscore,
        tsv_field(&d.numbered_seq),
    )
}

/// ANARCI-style per-position table: a comment header per domain, then one line
/// `<chain> <scheme-label> <aa>` for every scheme position (gaps shown as '-'), and
/// a `//` record terminator.
pub fn anarci_block(d: &DomainOut) -> String {
    let mut s = format!(
        "# {} | domain {} | {} | {} | {}..{}{}\n",
        d.query,
        d.domain,
        d.chain_type,
        d.scheme,
        d.start.map(|x| x.to_string()).unwrap_or_else(|| "?".into()),
        d.end,
        if d.complete { "" } else { " | partial" },
    );
    for (label, aa) in &d.all_positions {
        s.push_str(&format!("{} {} {}\n", d.chain_type, label, aa));
    }
    s.push_str("//\n");
    s
}
