//! Numbering schemes — a faithful port of ANARCI `schemes.py`.
//!
//! Shared engine (`smooth_insertions`, `number_regions`) + per-scheme functions.
//! Every function reproduces the reference behaviour byte-for-byte, including quirks.

use crate::constants::{alpha, blosum62, az_imgt, za_imgt, AZ26};
use crate::types::{assertion, CResult, Numbered, Res, State, StateType};

mod ig; // chothia / kabat / martin / wolfguy / aho live here

use StateType::{D, I, M};

// ---------------------------------------------------------------------------
// smooth_insertions — move HMMER framework-edge insertions into the adjacent CDR.
// ---------------------------------------------------------------------------

/// enforced_patterns[reg] (6 regions × 4 (state_id, type)).
const PATTERNS: [[(u8, StateType); 4]; 6] = [
    [(25, M), (26, M), (27, M), (28, I)],
    [(38, I), (38, M), (39, M), (40, M)],
    [(54, M), (55, M), (56, M), (57, I)],
    [(65, I), (65, M), (66, M), (67, M)],
    [(103, M), (104, M), (105, M), (106, I)],
    [(117, I), (117, M), (118, M), (119, M)],
];

pub fn smooth_insertions(state_vector: &[State]) -> Vec<State> {
    let mut buffer: Vec<State> = Vec::new();
    let mut sv: Vec<State> = Vec::new();
    let mut reg: i32 = -1; // only read when buffer is non-empty (always set first)

    for &st in state_vector {
        let id = st.id as i32;
        let buffered = if id < 23 {
            reg = -1;
            true
        } else if (25..28).contains(&id) {
            reg = 0;
            true
        } else if id > 37 && id <= 40 {
            reg = 1;
            true
        } else if (54..57).contains(&id) {
            reg = 2;
            true
        } else if id > 64 && id <= 67 {
            reg = 3;
            true
        } else if (103..106).contains(&id) {
            reg = 4;
            true
        } else if id > 116 && id <= 119 {
            reg = 5;
            true
        } else {
            false
        };

        if buffered {
            buffer.push(st);
            continue;
        }

        // Trigger state: flush the buffer (with pattern correction), then append.
        if !buffer.is_empty() {
            let nins = buffer.iter().filter(|s| s.typ == I).count();
            if nins > 0 {
                if reg == -1 {
                    flush_fw1(&buffer, nins, &mut sv);
                } else {
                    flush_cdr_edge(&buffer, reg as usize, &mut sv);
                }
            } else {
                sv.extend_from_slice(&buffer);
            }
            sv.push(st);
            buffer.clear();
        } else {
            sv.push(st);
        }
    }
    // NOTE: a trailing non-empty buffer is intentionally NOT flushed (matches ANARCI).
    sv
}

fn flush_fw1(buffer: &[State], nins: usize, sv: &mut Vec<State>) {
    let mut nt_dels = buffer[0].id as i32 - 1; // missing N-terminal states
    for s in buffer {
        if s.typ == D || s.si.is_none() {
            nt_dels += 1;
        } else {
            break; // first residue
        }
    }
    if nt_dels >= nins as i32 {
        let new_match: Vec<(u8, StateType)> =
            buffer.iter().filter(|s| s.typ == M).map(|s| (s.id, s.typ)).collect();
        let first = new_match[0].0 as i32;
        let nodel: Vec<State> = buffer.iter().copied().filter(|s| s.typ != D).collect();
        let add = nodel.len() as i32 - new_match.len() as i32;
        assert!(add >= 0, "Implementation logic error");
        let mut new_states: Vec<(u8, StateType)> = Vec::new();
        for x in (first - add)..first {
            new_states.push((x as u8, M));
        }
        new_states.extend(new_match);
        assert!(new_states.len() == nodel.len(), "Implementation logic error");
        for i in 0..nodel.len() {
            sv.push(State { id: new_states[i].0, typ: new_states[i].1, si: nodel[i].si });
        }
    } else {
        sv.extend_from_slice(buffer);
    }
}

fn flush_cdr_edge(buffer: &[State], reg: usize, sv: &mut Vec<State>) {
    let nodel: Vec<State> = buffer.iter().copied().filter(|s| s.typ != D).collect();
    let len = nodel.len();
    let pat = &PATTERNS[reg];
    let rep = len.saturating_sub(3);
    let new_states: Vec<(u8, StateType)> = if reg % 2 == 1 {
        // nterm fw: [pat[0]]*max(0,len-3) + pat[max(4-len,1):]
        let k = std::cmp::max(4i32 - len as i32, 1) as usize;
        let mut v: Vec<(u8, StateType)> = std::iter::repeat(pat[0]).take(rep).collect();
        v.extend_from_slice(&pat[k..]);
        v
    } else {
        // cterm fw: pat[:3] + [pat[2]]*max(0,len-3)
        let mut v: Vec<(u8, StateType)> = pat[..3].to_vec();
        v.extend(std::iter::repeat(pat[2]).take(rep));
        v
    };
    for i in 0..len {
        sv.push(State { id: new_states[i].0, typ: new_states[i].1, si: nodel[i].si });
    }
}

// ---------------------------------------------------------------------------
// number_regions — core region-splitting numbering engine.
// ---------------------------------------------------------------------------

/// Build the per-state region index array from a region string + char->index map.
pub(crate) fn region_map(region_string: &str, dict: &[(char, usize)]) -> Vec<usize> {
    region_string
        .chars()
        .map(|c| {
            dict.iter()
                .find(|(k, _)| *k == c)
                .map(|(_, v)| *v)
                .expect("region_string char must be in region_index_dict")
        })
        .collect()
}

type Regions = (Vec<Vec<Res>>, Option<usize>, Option<usize>);

#[allow(clippy::too_many_arguments)]
pub(crate) fn number_regions(
    sequence: &[u8],
    state_vector: &[State],
    state_string: &[u8],     // 'X' / 'I' per state, length 128
    region_of_state: &[usize], // region index per state, length 128
    rels: &mut [i32],
    n_regions: usize,
    exclude_deletions: &[usize],
) -> CResult<Regions> {
    let state_vector = smooth_insertions(state_vector);

    let mut regions: Vec<Vec<Res>> = vec![Vec::new(); n_regions];
    let mut insertion: i32 = -1;
    let mut previous_state_id: i32 = 1;
    let mut previous_state_type = D;
    let mut start_index: Option<usize> = None;
    let mut end_index: Option<usize> = None;
    let mut region: Option<usize> = None;

    for st in &state_vector {
        let sid = st.id as i32;
        let si_state = (st.id - 1) as usize;

        // Region selection — an insertion must NOT start a new region (BUG_FIX JD 9/4/15).
        if st.typ != I || region.is_none() {
            region = Some(region_of_state[si_state]);
        }
        let reg = region.unwrap();

        match st.typ {
            M => {
                if state_string[si_state] == b'I' {
                    if previous_state_type != D {
                        insertion += 1;
                    }
                    rels[reg] -= 1;
                } else {
                    insertion = -1;
                }
                regions[reg].push(Res::new(sid + rels[reg], alpha(insertion), sequence[st.si.unwrap()]));
                previous_state_id = sid;
                if start_index.is_none() {
                    start_index = st.si;
                }
                end_index = st.si;
                previous_state_type = M;
            }
            I => {
                insertion += 1;
                regions[reg].push(Res::new(
                    previous_state_id + rels[reg],
                    alpha(insertion),
                    sequence[st.si.unwrap()],
                ));
                if start_index.is_none() {
                    start_index = st.si;
                }
                end_index = st.si;
                previous_state_type = I;
            }
            D => {
                previous_state_type = D;
                if state_string[si_state] == b'I' {
                    rels[reg] -= 1;
                    continue; // skip overflow check (matches Python `continue`)
                }
                insertion = -1;
                previous_state_id = sid;
            }
        }

        if insertion >= 25 && exclude_deletions.contains(&reg) {
            insertion = 0;
        }
        if insertion >= 25 {
            return Err(assertion!("Too many insertions for numbering scheme to handle"));
        }
    }

    Ok((regions, start_index, end_index))
}

// ---------------------------------------------------------------------------
// CDR helpers
// ---------------------------------------------------------------------------

#[inline]
fn signed_index(i: isize, len: usize) -> usize {
    if i < 0 {
        (len as isize + i) as usize
    } else {
        i as usize
    }
}

/// `get_imgt_cdr`: symmetric CDR numbering (used by IMGT). Returns `len`-long
/// list of `Some((pos, ins))` or `None` (a gap the caller fills).
pub(crate) fn get_imgt_cdr(
    length: usize,
    maxlength: usize,
    start: i32,
    end: i32,
) -> Vec<Option<(i32, &'static str)>> {
    let n = length.max(maxlength);
    let mut ann: Vec<Option<(i32, &'static str)>> = vec![None; n];
    if length == 0 {
        return ann;
    }
    if length == 1 {
        ann[0] = Some((start, " "));
        return ann;
    }
    let mut front: isize = 0;
    let mut back: isize = -1;

    for i in 0..length.min(maxlength) {
        if i % 2 == 1 {
            ann[signed_index(back, n)] = Some((end + back as i32, " "));
            back -= 1;
        } else {
            ann[signed_index(front, n)] = Some((start + front as i32, " "));
            front += 1;
        }
    }

    let centre: Vec<usize> = (0..n).filter(|&i| ann[i].is_none()).collect();
    if centre.is_empty() {
        return ann;
    }
    let centre_left = ann[*centre.iter().min().unwrap() - 1].unwrap().0;
    let centre_right = ann[*centre.iter().max().unwrap() + 1].unwrap().0;

    let (frontfactor, backfactor): (isize, isize) = if maxlength % 2 == 0 {
        ((maxlength / 2) as isize, (maxlength / 2) as isize)
    } else {
        ((maxlength / 2 + 1) as isize, (maxlength / 2) as isize)
    };

    let extra = length.saturating_sub(maxlength);
    for i in 0..extra {
        if i % 2 == 0 {
            ann[signed_index(back, n)] = Some((centre_right, za52(back + backfactor)));
            back -= 1;
        } else {
            ann[signed_index(front, n)] = Some((centre_left, az52(front - frontfactor)));
            front += 1;
        }
    }
    ann
}

/// `az = alphabet[:-1]` (52 codes) with Python-style signed indexing.
#[inline]
fn az52(i: isize) -> &'static str {
    let idx = if i < 0 { (52 + i) as usize } else { i as usize };
    az_imgt(idx)
}
/// `za = az[::-1]` with Python-style signed indexing.
#[inline]
fn za52(i: isize) -> &'static str {
    let idx = if i < 0 { (52 + i) as usize } else { i as usize };
    za_imgt(idx)
}

/// `gap_missing`: fill skipped integer positions with gaps. Used by all schemes
/// except wolfguy. Input is the list of region lists.
pub(crate) fn gap_missing(numbering: &[Vec<Res>]) -> Vec<Res> {
    let mut num: Vec<Res> = vec![Res::gap(0)]; // sentinel ((0,' '),'-')
    for region in numbering {
        for &r in region {
            let last_num = num.last().unwrap().num;
            if r.num > last_num + 1 {
                for i in (last_num + 1)..r.num {
                    num.push(Res::gap(i));
                }
            }
            num.push(r);
        }
    }
    num[1..].to_vec()
}

/// `get_cdr3_annotations` for chothia/kabat (the branches actually used by the
/// schemes). Returns a sorted list of `(pos, ins)`. IMGT uses `get_imgt_cdr`.
pub(crate) fn get_cdr3_annotations(
    length: usize,
    scheme: &str,
    chain_type: &str,
) -> CResult<Vec<(i32, &'static str)>> {
    match (scheme, chain_type) {
        ("chothia", "heavy") | ("kabat", "heavy") => {
            let insertions = length.saturating_sub(10);
            if insertions >= 27 {
                return Err(assertion!("Too many insertions for numbering scheme to handle"));
            }
            let ordered: [(i32, &str); 10] = [
                (100, " "), (99, " "), (98, " "), (97, " "), (96, " "), (95, " "), (101, " "),
                (102, " "), (94, " "), (93, " "),
            ];
            let drop = 10usize.saturating_sub(length);
            let mut v: Vec<(i32, &str)> = ordered[drop..].to_vec();
            for a in &AZ26[..insertions] {
                v.push((100, a));
            }
            v.sort();
            Ok(v)
        }
        ("chothia", "light") | ("kabat", "light") => {
            let insertions = length.saturating_sub(9);
            if insertions >= 27 {
                return Err(assertion!("Too many insertions for numbering scheme to handle"));
            }
            let ordered: [(i32, &str); 9] = [
                (95, " "), (94, " "), (93, " "), (92, " "), (91, " "), (96, " "), (97, " "),
                (90, " "), (89, " "),
            ];
            let drop = 9usize.saturating_sub(length);
            let mut v: Vec<(i32, &str)> = ordered[drop..].to_vec();
            for a in &AZ26[..insertions] {
                v.push((95, a));
            }
            v.sort();
            Ok(v)
        }
        _ => Err(assertion!("Unimplemented scheme")),
    }
}

/// `_get_wolfguy_L1`: canonical-class CDRL1 numbering. Returns position list.
pub(crate) fn get_wolfguy_l1(seq: &[Res], length: usize) -> CResult<Vec<i32>> {
    // (name, consensus, positions) for each characterised length.
    let table: &[(usize, &[(&str, &[i32])])] = ig::WOLFGUY_L1;
    if let Some((_, canonicals)) = table.iter().find(|(l, _)| *l == length) {
        let mut best: Option<(&[i32], i32)> = None;
        for (consensus, positions) in canonicals.iter() {
            let cbytes = consensus.as_bytes();
            let mut sub_score = 0i32;
            for i in 0..length {
                let a = seq[i].aa.to_ascii_uppercase();
                let b = cbytes[i].to_ascii_uppercase();
                match blosum62(a, b) {
                    Some(s) => sub_score += s,
                    None => return Err(assertion!("Missing BLOSUM62 entry for wolfguy L1")),
                }
            }
            let better = match &best {
                Some((_, bs)) => sub_score > *bs, // strict: first canonical wins ties
                None => sub_score > -10000,
            };
            if better {
                best = Some((positions, sub_score));
            }
        }
        Ok(best.unwrap().0.to_vec())
    } else {
        // Symmetric numbering about the anchors.
        let mut ordered: Vec<i32> = Vec::new();
        for (p1, p2) in (551..575).zip((576..600).rev()) {
            ordered.push(p2);
            ordered.push(p1);
        }
        ordered.push(575);
        let mut v: Vec<i32> = ordered[..length.min(ordered.len())].to_vec();
        v.sort();
        Ok(v)
    }
}

// ---------------------------------------------------------------------------
// IMGT
// ---------------------------------------------------------------------------

pub fn number_imgt(state_vector: &[State], sequence: &[u8]) -> CResult<Numbered> {
    let state_string = b"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
    let region_string = "11111111111111111111111111222222222222333333333333333334444444444555555555555555555555555555555555555555666666666666677777777777";
    let dict = [('1', 0), ('2', 1), ('3', 2), ('4', 3), ('5', 4), ('6', 5), ('7', 6)];
    let region_of_state = region_map(region_string, &dict);
    let mut rels = [0i32; 7];
    let exclude_deletions = [1usize, 3, 5];

    let (regions, startindex, endindex) = number_regions(
        sequence,
        state_vector,
        state_string,
        &region_of_state,
        &mut rels,
        7,
        &exclude_deletions,
    )?;

    // CDR1: 27..39, maxlength 12.
    let cdr1 = regap_imgt_cdr(&regions[1], 12, 27, 39, 26);
    // CDR2: 56..66, maxlength 10.
    let cdr2 = regap_imgt_cdr(&regions[3], 10, 56, 66, 55);
    // CDR3: 105..118, maxlength 13.
    let cdr3seq: Vec<u8> = regions[5].iter().filter(|r| r.aa != b'-').map(|r| r.aa).collect();
    if cdr3seq.len() > 117 {
        return Ok(Numbered { residues: vec![], start: startindex, end: endindex });
    }
    let cdr3 = regap_imgt_cdr(&regions[5], 13, 105, 118, 104);

    let numbering = vec![
        regions[0].clone(),
        cdr1,
        regions[2].clone(),
        cdr2,
        regions[4].clone(),
        cdr3,
        regions[6].clone(),
    ];
    Ok(Numbered { residues: gap_missing(&numbering), start: startindex, end: endindex })
}

/// Shared CDR re-gapping used by IMGT (`get_imgt_cdr` + gap fill).
fn regap_imgt_cdr(region: &[Res], maxlength: usize, start: i32, end: i32, prev0: i32) -> Vec<Res> {
    let cdrseq: Vec<u8> = region.iter().filter(|r| r.aa != b'-').map(|r| r.aa).collect();
    let mut out: Vec<Res> = Vec::new();
    let mut si = 0usize;
    let mut prev_state = prev0;
    for ann in get_imgt_cdr(cdrseq.len(), maxlength, start, end) {
        match ann {
            None => {
                out.push(Res::gap(prev_state + 1));
                prev_state += 1;
            }
            Some((num, ins)) => {
                out.push(Res::new(num, ins, cdrseq[si]));
                prev_state = num;
                si += 1;
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Dispatcher
// ---------------------------------------------------------------------------

/// `number_sequence_from_alignment` — route to the right scheme/chain function.
/// Errors mirror ANARCI's `AssertionError(...)` messages exactly.
pub fn number_sequence_from_alignment(
    state_vector: &[State],
    sequence: &[u8],
    scheme: &str,
    chain_type: Option<&str>,
) -> CResult<Numbered> {
    let scheme = scheme.to_lowercase();
    let ct = chain_type.unwrap_or("");
    let is_kl = ct == "K" || ct == "L";
    match scheme.as_str() {
        "imgt" => number_imgt(state_vector, sequence),
        "chothia" => {
            if ct == "H" {
                ig::number_chothia_heavy(state_vector, sequence)
            } else if is_kl {
                ig::number_chothia_light(state_vector, sequence)
            } else {
                Err(unimplemented_err(&scheme, chain_type))
            }
        }
        "kabat" => {
            if ct == "H" {
                ig::number_kabat_heavy(state_vector, sequence)
            } else if is_kl {
                ig::number_kabat_light(state_vector, sequence)
            } else {
                Err(unimplemented_err(&scheme, chain_type))
            }
        }
        "martin" => {
            if ct == "H" {
                ig::number_martin_heavy(state_vector, sequence)
            } else if is_kl {
                ig::number_martin_light(state_vector, sequence)
            } else {
                Err(unimplemented_err(&scheme, chain_type))
            }
        }
        "aho" => ig::number_aho(state_vector, sequence, ct),
        "wolfguy" => {
            if ct == "H" {
                ig::number_wolfguy_heavy(state_vector, sequence)
            } else if is_kl {
                ig::number_wolfguy_light(state_vector, sequence)
            } else {
                Err(unimplemented_err(&scheme, chain_type))
            }
        }
        _ => Err(unimplemented_err(&scheme, chain_type)),
    }
}

fn unimplemented_err(scheme: &str, chain_type: Option<&str>) -> crate::types::CoreError {
    // Python: "Unimplemented numbering scheme %s for chain %s" % (scheme, chain_type)
    // chain_type may be None -> Python prints "None".
    let ct = match chain_type {
        Some(c) => c.to_string(),
        None => "None".to_string(),
    };
    assertion!("Unimplemented numbering scheme {} for chain {}", scheme, ct)
}
