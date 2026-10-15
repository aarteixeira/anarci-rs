//! Extract per-chain residues (sequence + Cα coordinate) from PDB files — robust to the
//! messy reality: multiple models, HETATM/MSE, alternate locations, waters, and missing
//! residues. One residue per Cα atom, in chain order; only model 1 is read.

use crate::structure::{parse_coord, ResCoord, StructChain};

pub fn three_to_one(res: &str) -> Option<u8> {
    Some(match res {
        "ALA" => b'A', "ARG" => b'R', "ASN" => b'N', "ASP" => b'D', "CYS" => b'C',
        "GLN" => b'Q', "GLU" => b'E', "GLY" => b'G', "HIS" => b'H', "ILE" => b'I',
        "LEU" => b'L', "LYS" => b'K', "MET" => b'M', "PHE" => b'F', "PRO" => b'P',
        "SER" => b'S', "THR" => b'T', "TRP" => b'W', "TYR" => b'Y', "VAL" => b'V',
        "MSE" => b'M', // selenomethionine — a common HETATM standing in for MET
        "SEC" => b'C', "PYL" => b'K',
        _ => return None,
    })
}

fn col(line: &[u8], a: usize, b: usize) -> &[u8] {
    let b = b.min(line.len());
    if a >= b { &[] } else { &line[a..b] }
}

/// Parse PDB into one [`StructChain`] per protein chain (id = `pdbid_chain`), with Cα coords.
pub fn chains_coords(pdb_bytes: &[u8]) -> Result<Vec<StructChain>, String> {
    let mut pdbid = String::from("pdb");
    let mut order: Vec<u8> = Vec::new();
    let mut res: std::collections::HashMap<u8, Vec<ResCoord>> = std::collections::HashMap::new();
    let mut seen: std::collections::HashSet<(u8, Vec<u8>)> = std::collections::HashSet::new();

    for raw in pdb_bytes.split(|&b| b == b'\n') {
        let line = raw.strip_suffix(b"\r").unwrap_or(raw);
        let rec = col(line, 0, 6);
        if rec.starts_with(b"HEADER") {
            let id = col(line, 62, 66).trim_ascii();
            if id.len() == 4 {
                pdbid = String::from_utf8_lossy(id).to_string();
            }
            continue;
        }
        if rec.starts_with(b"ENDMDL") {
            break; // model 1 only
        }
        if !(rec.starts_with(b"ATOM") || rec.starts_with(b"HETATM")) {
            continue;
        }
        if col(line, 12, 16).trim_ascii() != b"CA" {
            continue; // Cα only → one letter per residue
        }
        let alt = *col(line, 16, 17).first().unwrap_or(&b' ');
        if alt != b' ' && alt != b'A' {
            continue; // altloc: blank or primary
        }
        let resname = String::from_utf8_lossy(col(line, 17, 20).trim_ascii()).to_string();
        let Some(aa) = three_to_one(&resname) else { continue };
        let chain = *col(line, 21, 22).first().unwrap_or(&b' ');
        let key: Vec<u8> = col(line, 22, 27).to_vec(); // resSeq + iCode, deduped per chain
        if !seen.insert((chain, key)) {
            continue;
        }
        let ca = [parse_coord(col(line, 30, 38)), parse_coord(col(line, 38, 46)), parse_coord(col(line, 46, 54))];
        if !res.contains_key(&chain) {
            order.push(chain);
        }
        res.entry(chain).or_default().push(ResCoord { aa, ca });
    }

    Ok(order
        .into_iter()
        .map(|c| StructChain { id: format!("{pdbid}_{}", c as char), residues: res.remove(&c).unwrap_or_default() })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    // format one CA ATOM line at the correct fixed columns
    fn ca(serial: u32, res: &str, chain: char, resseq: i32, x: f32, y: f32, z: f32) -> String {
        format!("ATOM  {serial:>5}  CA  {res} {chain}{resseq:>4}    {x:>8.3}{y:>8.3}{z:>8.3}  1.00  0.00           C")
    }

    #[test]
    fn extracts_ca_sequence_and_coords() {
        let mut pdb = String::from("HEADER    IMMUNOGLOBULIN                          01-JAN-00   1ABC\n");
        pdb += "ATOM      1  N   GLU H   1       0.000   0.000   0.000  1.00  0.00           N\n"; // non-CA, ignored
        pdb += &(ca(2, "GLU", 'H', 1, 1.0, 2.0, 3.0) + "\n");
        pdb += &(ca(3, "VAL", 'H', 2, 4.0, 5.0, 6.0) + "\n");
        pdb += "ATOM      9  CA  HOH H 101       9.000   9.000   9.000  1.00  0.00           O\n"; // water, skipped
        pdb += "ENDMDL\n";
        let cs = chains_coords(pdb.as_bytes()).unwrap();
        assert_eq!(cs.len(), 1);
        assert_eq!(cs[0].id, "1ABC_H");
        assert_eq!(cs[0].seq(), b"EV");
        assert_eq!(cs[0].residues[0].ca, [1.0, 2.0, 3.0]);
        assert_eq!(cs[0].residues[1].ca, [4.0, 5.0, 6.0]);
    }
}
