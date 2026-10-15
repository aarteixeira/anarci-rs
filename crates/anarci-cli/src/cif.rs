//! Extract first-model protein C-alpha coordinates from an mmCIF atom-site loop.
//! Quoted fields and wrapped rows are supported. Semicolon text fields are tokenized.
use crate::pdb::three_to_one;
use crate::structure::{parse_coord, ResCoord, StructChain};

struct Token {
    value: String,
    quoted: bool,
}
fn tokens(text: &str) -> Result<Vec<Token>, String> {
    let b = text.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        if b[i].is_ascii_whitespace() {
            i += 1;
            continue;
        }
        if b[i] == b'#' {
            while i < b.len() && b[i] != b'\n' {
                i += 1;
            }
            continue;
        }
        let start = i;
        if b[i] == b';' && (i == 0 || b[i - 1] == b'\n') {
            i += 1;
            let content = i;
            while i < b.len() && !(b[i] == b';' && b[i - 1] == b'\n') {
                i += 1;
            }
            if i == b.len() {
                return Err("mmCIF: unterminated semicolon text field".into());
            }
            out.push(Token {
                value: text[content..i].to_string(),
                quoted: true,
            });
            i += 1;
        } else if b[i] == b'\'' || b[i] == b'"' {
            let quote = b[i];
            i += 1;
            let content = i;
            while i < b.len()
                && !(b[i] == quote && (i + 1 == b.len() || b[i + 1].is_ascii_whitespace()))
            {
                i += 1;
            }
            if i == b.len() {
                return Err("mmCIF: unterminated quoted field".into());
            }
            out.push(Token {
                value: text[content..i].to_string(),
                quoted: true,
            });
            i += 1;
        } else {
            while i < b.len() && !b[i].is_ascii_whitespace() {
                i += 1;
            }
            out.push(Token {
                value: text[start..i].to_string(),
                quoted: false,
            });
        }
    }
    Ok(out)
}
fn control(t: &Token) -> bool {
    !t.quoted
        && (t.value == "loop_"
            || t.value == "stop_"
            || t.value.starts_with('_')
            || t.value.starts_with("data_")
            || t.value.starts_with("save_"))
}

pub fn chains_coords(cif: &[u8]) -> Result<Vec<StructChain>, String> {
    let text = std::str::from_utf8(cif).map_err(|_| "mmCIF: input must be UTF-8")?;
    let tokens = tokens(text)?;
    let mut i = 0;
    let mut cols = Vec::new();
    let mut rows = None;
    while i < tokens.len() {
        if !tokens[i].quoted && tokens[i].value == "loop_" {
            i += 1;
            cols.clear();
            while i < tokens.len() && !tokens[i].quoted && tokens[i].value.starts_with('_') {
                cols.push(tokens[i].value.clone());
                i += 1;
            }
            if cols.iter().any(|s| s.starts_with("_atom_site.")) {
                if cols.iter().any(|s| !s.starts_with("_atom_site.")) {
                    return Err("mmCIF: mixed categories in atom-site loop".into());
                }
                let start = i;
                while i < tokens.len() && !control(&tokens[i]) {
                    i += 1;
                }
                rows = Some(&tokens[start..i]);
                break;
            }
        } else {
            i += 1;
        }
    }
    let rows = rows.ok_or("mmCIF: no _atom_site loop found")?;
    if rows.len() % cols.len() != 0 {
        return Err("mmCIF: incomplete atom-site row".into());
    }
    let idx = |name: &str| cols.iter().position(|c| c == &format!("_atom_site.{name}"));
    let c_group = idx("group_PDB");
    let c_atom = idx("label_atom_id").ok_or("mmCIF: label_atom_id missing")?;
    let c_comp = idx("label_comp_id").ok_or("mmCIF: label_comp_id missing")?;
    let c_label_chain = idx("label_asym_id");
    let c_chain = idx("auth_asym_id")
        .or_else(|| idx("label_asym_id"))
        .ok_or("mmCIF: no chain id column")?;
    let c_label_seq = idx("label_seq_id");
    let c_seq = idx("auth_seq_id")
        .or_else(|| idx("label_seq_id"))
        .ok_or("mmCIF: residue sequence id column missing")?;
    let c_ins = idx("pdbx_PDB_ins_code");
    let c_alt = idx("label_alt_id");
    let c_model = idx("pdbx_PDB_model_num");
    let (c_x, c_y, c_z) = (idx("Cartn_x"), idx("Cartn_y"), idx("Cartn_z"));
    let mut order = Vec::new();
    let mut res: std::collections::HashMap<String, Vec<ResCoord>> =
        std::collections::HashMap::new();
    let mut seen = std::collections::HashSet::new();
    let mut first_model: Option<&str> = None;
    for row in rows.chunks(cols.len()) {
        let field = |c: usize| row[c].value.as_str();
        if let Some(g) = c_group {
            if !matches!(field(g), "ATOM" | "HETATM") {
                continue;
            }
        }
        if field(c_atom) != "CA" {
            continue;
        }
        if let Some(m) = c_model {
            match first_model {
                None => first_model = Some(field(m)),
                Some(first) if first != field(m) => continue,
                _ => {}
            }
        }
        if let Some(a) = c_alt {
            if !matches!(field(a), "." | "?" | "A") {
                continue;
            }
        }
        let Some(aa) = three_to_one(field(c_comp)) else {
            continue;
        };
        let preferred = |primary, fallback: Option<usize>| {
            let value = field(primary);
            if matches!(value, "?" | ".") {
                fallback.map(field).unwrap_or(value)
            } else {
                value
            }
        };
        let chain = preferred(c_chain, c_label_chain).to_string();
        let seq_id = preferred(c_seq, c_label_seq);
        if matches!(seq_id, "?" | ".") {
            return Err("mmCIF: protein C-alpha residue lacks a sequence identifier".into());
        }
        if !seen.insert((
            chain.clone(),
            seq_id.to_string(),
            c_ins.map(field).unwrap_or("?").to_string(),
        )) {
            continue;
        }
        let cf = |c: Option<usize>| {
            c.map(|c| parse_coord(field(c).as_bytes()))
                .unwrap_or(f32::NAN)
        };
        if !res.contains_key(&chain) {
            order.push(chain.clone());
        }
        res.entry(chain).or_default().push(ResCoord {
            aa,
            ca: [cf(c_x), cf(c_y), cf(c_z)],
        });
    }
    Ok(order
        .into_iter()
        .map(|c| StructChain {
            id: format!("cif_{c}"),
            residues: res.remove(&c).unwrap_or_default(),
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extracts_ca_from_atom_site_loop() {
        let cif = "\
data_test
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_alt_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.pdbx_PDB_model_num
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM N   GLU . H 1 ? 1 0 0 0
ATOM CA  GLU . H 1 ? 1 1 2 3
ATOM CB  GLU . H 1 ? 1 0 0 0
ATOM CA  VAL A H 2 ? 1 4 5 6
ATOM CA  VAL B H 2 ? 1 9 9 9
HETATM CA  HOH . H 101 ? 1 0 0 0
ATOM CA  GLY . H 3 ? 2 0 0 0
#
";
        let cs = chains_coords(cif.as_bytes()).unwrap();
        assert_eq!(cs.len(), 1);
        assert_eq!(cs[0].id, "cif_H");
        // GLU + VAL(altloc A) only; altloc B dropped, water skipped, model-2 GLY excluded
        assert_eq!(cs[0].seq(), b"EV");
        assert_eq!(cs[0].residues[0].ca, [1.0, 2.0, 3.0]);
        assert_eq!(cs[0].residues[1].ca, [4.0, 5.0, 6.0]);
    }
    #[test]
    fn quoted_chains_wrapped_rows_and_malformed_fields() {
        let header = "data_synthetic\nloop_\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";
        let cs = chains_coords(
            format!("{header}'CA' ALA 'chain one' 1\n1 2 3\nCA GLY 'chain one' 2 4 5 6\n#\n")
                .as_bytes(),
        )
        .unwrap();
        assert_eq!(cs[0].id, "cif_chain one");
        assert_eq!(cs[0].seq(), b"AG");
        assert_eq!(cs[0].residues[1].ca, [4., 5., 6.]);
        assert!(chains_coords(format!("{header}CA ALA A 1 2\n").as_bytes()).is_err());
        assert!(chains_coords(format!("{header}'CA ALA A 1 1 2 3\n").as_bytes()).is_err());
    }
}
