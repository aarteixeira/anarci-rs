//! Streaming FASTA / FASTQ readers (low memory for NGS-scale inputs) and
//! nucleotide/protein detection.

use std::io::BufRead;

/// One input sequence record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Record {
    pub id: String,
    pub desc: String,
    pub seq: Vec<u8>,
}

impl Record {
    /// `id` plus any description, as one label (what ANARCI calls the query name).
    pub fn name(&self) -> String {
        if self.desc.is_empty() {
            self.id.clone()
        } else {
            format!("{} {}", self.id, self.desc)
        }
    }
}

/// FASTA or FASTQ, decided by the first non-blank byte (`>` vs `@`).
pub fn read_records<R: BufRead>(mut r: R) -> Result<Vec<Record>, String> {
    loop {
        let buf = r.fill_buf().map_err(|e| e.to_string())?;
        if buf.is_empty() {
            return Ok(Vec::new());
        }
        let skip = buf.iter().take_while(|b| b.is_ascii_whitespace()).count();
        if skip == 0 {
            break;
        }
        r.consume(skip);
    }
    match r.fill_buf().map_err(|e| e.to_string())?[0] {
        b'>' => read_fasta(r),
        b'@' => read_fastq(r),
        c => Err(format!(
            "unrecognized input: first non-blank byte {:?} is neither FASTA '>' nor FASTQ '@'",
            c as char
        )),
    }
}

/// Streaming FASTA: `>id desc` then wrapped sequence lines until the next `>`.
pub fn read_fasta<R: BufRead>(r: R) -> Result<Vec<Record>, String> {
    let mut out = Vec::new();
    let mut cur: Option<Record> = None;
    for line in r.lines() {
        let line = line.map_err(|e| e.to_string())?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(rec) = cur.take() {
                out.push(rec);
            }
            let (id, desc) = split_header(header);
            if id.is_empty() {
                return Err("sequence identifier must not be empty".into());
            }
            cur = Some(Record {
                id,
                desc,
                seq: Vec::new(),
            });
        } else if let Some(rec) = cur.as_mut() {
            rec.seq
                .extend(line.bytes().filter(|b| !b.is_ascii_whitespace()));
        } else {
            return Err("FASTA sequence occurs before its header".into());
        }
    }
    if let Some(rec) = cur.take() {
        out.push(rec);
    }
    Ok(out)
}

/// Streaming FASTQ: 4-line records (`@id`, seq, `+`, qual). Quality is discarded.
pub fn read_fastq<R: BufRead>(r: R) -> Result<Vec<Record>, String> {
    let mut out = Vec::new();
    let mut lines = r.lines();
    while let Some(h) = lines.next() {
        let h = h.map_err(|e| e.to_string())?;
        if h.trim().is_empty() {
            continue;
        }
        let header = h
            .strip_prefix('@')
            .ok_or_else(|| format!("FASTQ record must start with '@', got {h:?}"))?;
        let (id, desc) = split_header(header);
        let seq = lines
            .next()
            .ok_or("FASTQ truncated: missing sequence line")?
            .map_err(|e| e.to_string())?;
        let plus = lines
            .next()
            .ok_or("FASTQ truncated: missing '+' line")?
            .map_err(|e| e.to_string())?;
        let qual = lines
            .next()
            .ok_or("FASTQ truncated: missing quality line")?
            .map_err(|e| e.to_string())?;
        if id.is_empty() {
            return Err("FASTQ identifier must not be empty".into());
        }
        let repeated = plus
            .strip_prefix('+')
            .ok_or("FASTQ separator must start with '+'")?;
        if !repeated.is_empty() && split_header(repeated).0 != id {
            return Err(format!(
                "FASTQ separator identifier differs from header: {id}"
            ));
        }
        if seq.len() != qual.len() {
            return Err(format!("FASTQ sequence/quality length mismatch for {id}"));
        }
        if seq.is_empty() || seq.bytes().any(|b| b.is_ascii_whitespace()) {
            return Err(format!(
                "FASTQ sequence must be a non-empty line without whitespace: {id}"
            ));
        }
        if qual.bytes().any(|b| !(33..=126).contains(&b)) {
            return Err(format!(
                "FASTQ quality contains invalid ASCII characters: {id}"
            ));
        }
        out.push(Record {
            id,
            desc,
            seq: seq.into_bytes(),
        });
    }
    Ok(out)
}

fn split_header(h: &str) -> (String, String) {
    let h = h.trim();
    match h.split_once(char::is_whitespace) {
        Some((id, desc)) => (id.to_string(), desc.trim().to_string()),
        None => (h.to_string(), String::new()),
    }
}

/// Heuristic: a sequence is nucleotide if ≥90% of its letters are in {A,C,G,T,U,N}.
/// (Ambiguity codes beyond N are rare in reads and fall through to protein, which
/// is the safe default — a mis-called protein just fails to number.)
pub fn is_nucleotide(seq: &[u8]) -> bool {
    let mut acgt = 0usize;
    let mut letters = 0usize;
    for &b in seq {
        if b.is_ascii_alphabetic() {
            letters += 1;
            if matches!(
                b.to_ascii_uppercase(),
                b'A' | b'C' | b'G' | b'T' | b'U' | b'N'
            ) {
                acgt += 1;
            }
        }
    }
    letters > 0 && (acgt as f64 / letters as f64) >= 0.9
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn blank_prefix_and_wrapped_fasta() {
        let records = read_records(&b"\n \n>r description\nAC D\nEF\n"[..]).unwrap();
        assert_eq!(records[0].name(), "r description");
        assert_eq!(records[0].seq, b"ACDEF");
        assert!(read_records(&b" \n\t"[..]).unwrap().is_empty());
    }
    #[test]
    fn validates_fastq() {
        assert_eq!(
            read_records(&b"\n@r\nACG\n+r\n!!!\n"[..]).unwrap()[0].seq,
            b"ACG"
        );
        for invalid in [
            "@r\nAC\n-\n!!\n",
            "@r\nAC\n+\n!\n",
            "@r\nAC\n+different\n!!\n",
            "@r\nAC\n+\n",
            "@\nAC\n+\n!!\n",
        ] {
            assert!(read_records(invalid.as_bytes()).is_err(), "{invalid}");
        }
        assert!(read_fasta(&b"AC\n>r\nDE\n"[..]).is_err());
        assert!(read_records(&b">\nAC\n"[..]).is_err());
    }
}
