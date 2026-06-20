//! Chothia / Kabat / Martin / Wolfguy / Aho scheme functions.
//!
//! Faithful port of ANARCI `schemes.py` (functions verified byte-correct against
//! the golden oracle). Each builds the per-scheme inputs and calls the shared
//! `number_regions` engine in `super`, then re-gaps the CDR/FW regions exactly as
//! the reference does.

use crate::constants::alpha;
use crate::types::{assertion, CResult, Numbered, Res, State};

/// Wolfguy CDRL1 canonical-class table: (length, [(consensus, positions)]).
/// Names from schemes.py are dropped (only consensus + positions are used).
pub(crate) static WOLFGUY_L1: &[(usize, &[(&str, &[i32])])] = &[
    (9, &[("XXXXXXXXX", &[551, 552, 554, 556, 563, 572, 597, 598, 599])]),
    (10, &[("XXXXXXXXXX", &[551, 552, 553, 556, 561, 562, 571, 597, 598, 599])]),
    (
        11,
        &[
            ("RASQDISSYLA", &[551, 552, 553, 556, 561, 562, 571, 596, 597, 598, 599]),
            ("GGNNIGSKSVH", &[551, 552, 554, 556, 561, 562, 571, 572, 597, 598, 599]),
            ("SGDQLPKKYAY", &[551, 552, 554, 556, 561, 562, 571, 572, 597, 598, 599]),
        ],
    ),
    (
        12,
        &[
            ("TLSSQHSTYTIE", &[551, 552, 553, 554, 555, 556, 561, 563, 572, 597, 598, 599]),
            ("TASSSVSSSYLH", &[551, 552, 553, 556, 561, 562, 571, 595, 596, 597, 598, 599]),
            ("RASQSVxNNYLA", &[551, 552, 553, 556, 561, 562, 571, 581, 596, 597, 598, 599]),
            ("rSShSIrSrrVh", &[551, 552, 553, 556, 561, 562, 571, 581, 596, 597, 598, 599]),
        ],
    ),
    (
        13,
        &[
            ("SGSSSNIGNNYVS", &[551, 552, 554, 555, 556, 557, 561, 562, 571, 572, 597, 598, 599]),
            ("TRSSGSLANYYVQ", &[551, 552, 553, 554, 556, 561, 562, 563, 571, 572, 597, 598, 599]),
        ],
    ),
    (
        14,
        &[
            ("RSSTGAVTTSNYAN", &[551, 552, 553, 554, 555, 561, 562, 563, 564, 571, 572, 597, 598, 599]),
            ("TGTSSDVGGYNYVS", &[551, 552, 554, 555, 556, 557, 561, 562, 571, 572, 596, 597, 598, 599]),
        ],
    ),
    (15, &[("XXXXXXXXXXXXXXX", &[551, 552, 553, 556, 561, 562, 563, 581, 582, 594, 595, 596, 597, 598, 599])]),
    (16, &[("XXXXXXXXXXXXXXXX", &[551, 552, 553, 556, 561, 562, 563, 581, 582, 583, 594, 595, 596, 597, 598, 599])]),
    (17, &[("XXXXXXXXXXXXXXXXX", &[551, 552, 553, 556, 561, 562, 563, 581, 582, 583, 584, 594, 595, 596, 597, 598, 599])]),
];

// ---------------------------------------------------------------------------
// Helpers shared by the ports.
// ---------------------------------------------------------------------------

/// Zip a list of `(pos, ins)` annotations with the residues of a region, producing
/// `Res` for each. Mirrors Python `[(annotations[i], region[i][1]) for i in range(length)]`.
fn apply_annotations(annotations: &[(i32, &'static str)], region: &[Res]) -> Vec<Res> {
    (0..region.len())
        .map(|i| Res::new(annotations[i].0, annotations[i].1, region[i].aa))
        .collect()
}

// ---------------------------------------------------------------------------
// Chothia — heavy
// ---------------------------------------------------------------------------

pub(crate) fn number_chothia_heavy(state_vector: &[State], sequence: &[u8]) -> CResult<Numbered> {
    let state_string = b"XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX";
    let region_string = "11111111112222222222222333333333333333444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888";
    let dict = [('1', 0), ('2', 1), ('3', 2), ('4', 3), ('5', 4), ('6', 5), ('7', 6), ('8', 7)];
    let region_of_state = super::region_map(region_string, &dict);
    let mut rels = [0i32, -1, -1, -5, -5, -8, -12, -15];
    let exclude_deletions = [0usize, 2, 4, 6];

    let (regions, startindex, endindex) = super::number_regions(
        sequence, state_vector, state_string, &region_of_state, &mut rels, 8, &exclude_deletions,
    )?;

    let mut _numbering: Vec<Vec<Res>> = vec![Vec::new(); 8];
    _numbering[1] = regions[1].clone();
    _numbering[3] = regions[3].clone();
    _numbering[5] = regions[5].clone();
    _numbering[7] = regions[7].clone();

    // Region 1 (index 0): insertions placed at Chothia position 6.
    let insertions = regions[0].iter().filter(|r| r.ins != " ").count();
    if insertions > 0 {
        let start = regions[0][0].num;
        let length = regions[0].len();
        let mut ann: Vec<(i32, &str)> = Vec::new();
        for p in start..7 {
            ann.push((p, " "));
        }
        for i in 0..insertions {
            ann.push((6, alpha(i as i32)));
        }
        ann.push((7, " "));
        ann.push((8, " "));
        ann.push((9, " "));
        _numbering[0] = apply_annotations(&ann[..length], &regions[0]);
    } else {
        _numbering[0] = regions[0].clone();
    }

    // CDR1 (index 2): insertions onto 31.
    _numbering[2] = cdrh1_chothia_martin(&regions[2]);

    // CDR2 (index 4): insertions onto 52.
    _numbering[4] = cdrh2_50_57(&regions[4]);

    // CDR3 (index 6): insertions onto 100.
    if regions[6].len() > 36 {
        return Ok(Numbered { residues: vec![], start: startindex, end: endindex });
    }
    let ann = super::get_cdr3_annotations(regions[6].len(), "chothia", "heavy")?;
    _numbering[6] = apply_annotations(&ann, &regions[6]);

    Ok(Numbered { residues: super::gap_missing(&_numbering), start: startindex, end: endindex })
}

/// CDRH1 23..33 gapping (Chothia / Martin heavy, insertions on 31).
fn cdrh1_chothia_martin(region: &[Res]) -> Vec<Res> {
    let length = region.len();
    let insertions = length.saturating_sub(11);
    let mut ann: Vec<(i32, &str)> = Vec::new();
    if insertions > 0 {
        for p in 23..32 {
            ann.push((p, " "));
        }
        for i in 0..insertions {
            ann.push((31, alpha(i as i32)));
        }
        ann.push((32, " "));
        ann.push((33, " "));
    } else {
        // [(23..32)][:length-2] + [(32),(33)][:length]
        let base: Vec<(i32, &str)> = (23..32).map(|p| (p, " ")).collect();
        let take = (length as i32 - 2).max(0) as usize;
        ann.extend_from_slice(&base[..take.min(base.len())]);
        let tail = [(32i32, " "), (33, " ")];
        ann.extend_from_slice(&tail[..length.min(2)]);
    }
    apply_annotations(&ann, region)
}

/// CDRH2 50..57 gapping shared by chothia/kabat/martin heavy (insertions on 52).
fn cdrh2_50_57(region: &[Res]) -> Vec<Res> {
    let length = region.len();
    let insertions = length.saturating_sub(8);
    let mut ann: Vec<(i32, &str)> = Vec::new();
    // [(50),(51),(52)][:max(0,length-5)]
    let head = [(50i32, " "), (51, " "), (52, " ")];
    let head_take = (length as i32 - 5).max(0) as usize;
    ann.extend_from_slice(&head[..head_take.min(3)]);
    for i in 0..insertions {
        ann.push((52, alpha(i as i32)));
    }
    // [(53),(54),(55),(56),(57)][abs(min(0,length-5)):]
    let tail = [(53i32, " "), (54, " "), (55, " "), (56, " "), (57, " ")];
    let drop = (length as i32 - 5).min(0).unsigned_abs() as usize;
    ann.extend_from_slice(&tail[drop.min(5)..]);
    apply_annotations(&ann, region)
}

// ---------------------------------------------------------------------------
// Chothia — light
// ---------------------------------------------------------------------------

pub(crate) fn number_chothia_light(state_vector: &[State], sequence: &[u8]) -> CResult<Numbered> {
    let state_string = b"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX";
    let region_string = "11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666666677777777777";
    let dict = [('1', 0), ('2', 1), ('3', 2), ('4', 3), ('5', 4), ('6', 5), ('7', 6)];
    let region_of_state = super::region_map(region_string, &dict);
    let mut rels = [0i32, 0, -6, -6, -13, -16, -20];
    let exclude_deletions = [1usize, 3, 4, 5];

    let (regions, startindex, endindex) = super::number_regions(
        sequence, state_vector, state_string, &region_of_state, &mut rels, 7, &exclude_deletions,
    )?;

    let mut _numbering: Vec<Vec<Res>> = vec![Vec::new(); 7];
    _numbering[0] = regions[0].clone();
    _numbering[2] = regions[2].clone();
    _numbering[4] = regions[4].clone();
    _numbering[6] = regions[6].clone();

    // CDR1 (index 1): insertions onto 30.
    {
        let length = regions[1].len();
        let insertions = length.saturating_sub(11);
        let head = [(24i32, " "), (25, " "), (26, " "), (27, " "), (28, " "), (29, " "), (30, " ")];
        let mut ann: Vec<(i32, &str)> = head[..length.min(7)].to_vec();
        for i in 0..insertions {
            ann.push((30, alpha(i as i32)));
        }
        let tail = [(31i32, " "), (32, " "), (33, " "), (34, " ")];
        let drop = (length as i32 - 11).min(0).unsigned_abs() as usize;
        ann.extend_from_slice(&tail[drop.min(4)..]);
        _numbering[1] = apply_annotations(&ann, &regions[1]);
    }

    // CDR2 (index 3): insertions onto 52.
    _numbering[3] = cdrl2_chothia(&regions[3]);

    // FW3 (index 4): insertions on 68; first deletion on 68 at length 33; else alignment.
    {
        let length = regions[4].len();
        let insertions = length.saturating_sub(34);
        if insertions > 0 {
            let mut ann: Vec<(i32, &str)> = (55..69).map(|p| (p, " ")).collect();
            for i in 0..insertions {
                ann.push((68, alpha(i as i32)));
            }
            ann.extend((69..89).map(|p| (p, " ")));
            _numbering[4] = apply_annotations(&ann, &regions[4]);
        } else if length == 33 {
            let mut ann: Vec<(i32, &str)> = (55..68).map(|p| (p, " ")).collect();
            ann.extend((69..89).map(|p| (p, " ")));
            _numbering[4] = apply_annotations(&ann, &regions[4]);
        } else {
            _numbering[4] = regions[4].clone();
        }
    }

    // CDR3 (index 5): insertions onto 95.
    if regions[5].len() > 35 {
        return Ok(Numbered { residues: vec![], start: startindex, end: endindex });
    }
    let ann = super::get_cdr3_annotations(regions[5].len(), "chothia", "light")?;
    _numbering[5] = apply_annotations(&ann, &regions[5]);

    Ok(Numbered { residues: super::gap_missing(&_numbering), start: startindex, end: endindex })
}

/// CDRL2 51..54 gapping shared by chothia/kabat/martin light (insertions on 52).
fn cdrl2_chothia(region: &[Res]) -> Vec<Res> {
    let length = region.len();
    let insertions = length.saturating_sub(4);
    if insertions > 0 {
        let mut ann: Vec<(i32, &str)> = vec![(51, " "), (52, " ")];
        for i in 0..insertions {
            ann.push((52, alpha(i as i32)));
        }
        ann.push((53, " "));
        ann.push((54, " "));
        apply_annotations(&ann, region)
    } else {
        // alignment placement
        region.to_vec()
    }
}

// ---------------------------------------------------------------------------
// Kabat — heavy
// ---------------------------------------------------------------------------

pub(crate) fn number_kabat_heavy(state_vector: &[State], sequence: &[u8]) -> CResult<Numbered> {
    let state_string = b"XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX";
    let region_string = "11111111112222222222222333333333333333334444444444444455555555555666666666666666666666666666666666666666777777777777788888888888";
    let dict = [('1', 0), ('2', 1), ('3', 2), ('4', 3), ('5', 4), ('6', 5), ('7', 6), ('8', 7)];
    let region_of_state = super::region_map(region_string, &dict);
    let mut rels = [0i32, -1, -1, -5, -5, -8, -12, -15];
    let exclude_deletions = [2usize, 4, 6];

    let (regions, startindex, endindex) = super::number_regions(
        sequence, state_vector, state_string, &region_of_state, &mut rels, 8, &exclude_deletions,
    )?;

    let mut _numbering: Vec<Vec<Res>> = vec![Vec::new(); 8];
    _numbering[1] = regions[1].clone();
    _numbering[3] = regions[3].clone();
    _numbering[5] = regions[5].clone();
    _numbering[7] = regions[7].clone();

    // Region 1 (index 0): insertions placed at Kabat position 6.
    let insertions = regions[0].iter().filter(|r| r.ins != " ").count();
    if insertions > 0 {
        let start = regions[0][0].num;
        let length = regions[0].len();
        let mut ann: Vec<(i32, &str)> = Vec::new();
        for p in start..7 {
            ann.push((p, " "));
        }
        for i in 0..insertions {
            ann.push((6, alpha(i as i32)));
        }
        ann.push((7, " "));
        ann.push((8, " "));
        ann.push((9, " "));
        _numbering[0] = apply_annotations(&ann[..length], &regions[0]);
    } else {
        _numbering[0] = regions[0].clone();
    }

    // CDR1 (index 2): insertions onto 35, delete from 35 backwards.
    {
        let length = regions[2].len();
        let insertions = length.saturating_sub(13);
        let base: Vec<(i32, &str)> = (23..36).map(|p| (p, " ")).collect();
        let mut ann: Vec<(i32, &str)> = base[..length.min(base.len())].to_vec();
        for i in 0..insertions {
            ann.push((35, alpha(i as i32)));
        }
        _numbering[2] = apply_annotations(&ann, &regions[2]);
    }

    // CDR2 (index 4): insertions onto 52.
    _numbering[4] = cdrh2_50_57(&regions[4]);

    // CDR3 (index 6): insertions onto 100.
    if regions[6].len() > 36 {
        return Ok(Numbered { residues: vec![], start: startindex, end: endindex });
    }
    let ann = super::get_cdr3_annotations(regions[6].len(), "kabat", "heavy")?;
    _numbering[6] = apply_annotations(&ann, &regions[6]);

    Ok(Numbered { residues: super::gap_missing(&_numbering), start: startindex, end: endindex })
}

// ---------------------------------------------------------------------------
// Kabat — light
// ---------------------------------------------------------------------------

pub(crate) fn number_kabat_light(state_vector: &[State], sequence: &[u8]) -> CResult<Numbered> {
    let state_string = b"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX";
    let region_string = "11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666666677777777777";
    let dict = [('1', 0), ('2', 1), ('3', 2), ('4', 3), ('5', 4), ('6', 5), ('7', 6)];
    let region_of_state = super::region_map(region_string, &dict);
    let mut rels = [0i32, 0, -6, -6, -13, -16, -20];
    let exclude_deletions = [1usize, 3, 5];

    let (regions, startindex, endindex) = super::number_regions(
        sequence, state_vector, state_string, &region_of_state, &mut rels, 7, &exclude_deletions,
    )?;

    let mut _numbering: Vec<Vec<Res>> = vec![Vec::new(); 7];
    _numbering[0] = regions[0].clone();
    _numbering[2] = regions[2].clone();
    _numbering[4] = regions[4].clone();
    _numbering[6] = regions[6].clone();

    // CDR1 (index 1): insertions onto 27, delete forward from 28.
    {
        let length = regions[1].len();
        let insertions = length.saturating_sub(11);
        let head = [(24i32, " "), (25, " "), (26, " "), (27, " ")];
        let mut ann: Vec<(i32, &str)> = head[..length.min(4)].to_vec();
        for i in 0..insertions {
            ann.push((27, alpha(i as i32)));
        }
        let tail = [(28i32, " "), (29, " "), (30, " "), (31, " "), (32, " "), (33, " "), (34, " ")];
        let drop = (length as i32 - 11).min(0).unsigned_abs() as usize;
        ann.extend_from_slice(&tail[drop.min(7)..]);
        _numbering[1] = apply_annotations(&ann, &regions[1]);
    }

    // CDR2 (index 3): insertions onto 52.
    _numbering[3] = cdrl2_chothia(&regions[3]);

    // FW3: all insertions placed by alignment (region 4 left unchanged).
    // CDR3 (index 5): insertions onto 95.
    if regions[5].len() > 35 {
        return Ok(Numbered { residues: vec![], start: startindex, end: endindex });
    }
    let ann = super::get_cdr3_annotations(regions[5].len(), "kabat", "light")?;
    _numbering[5] = apply_annotations(&ann, &regions[5]);

    Ok(Numbered { residues: super::gap_missing(&_numbering), start: startindex, end: endindex })
}

// ---------------------------------------------------------------------------
// Martin (extended Chothia) — heavy
// ---------------------------------------------------------------------------

pub(crate) fn number_martin_heavy(state_vector: &[State], sequence: &[u8]) -> CResult<Numbered> {
    let state_string = b"XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXIIIXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX";
    let region_string = "11111111112222222222222333333333333333444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888";
    let dict = [('1', 0), ('2', 1), ('3', 2), ('4', 3), ('5', 4), ('6', 5), ('7', 6), ('8', 7)];
    let region_of_state = super::region_map(region_string, &dict);
    let mut rels = [0i32, -1, -1, -5, -5, -8, -12, -15];
    let exclude_deletions = [2usize, 4, 5, 6];

    let (regions, startindex, endindex) = super::number_regions(
        sequence, state_vector, state_string, &region_of_state, &mut rels, 8, &exclude_deletions,
    )?;

    let mut _numbering: Vec<Vec<Res>> = vec![Vec::new(); 8];
    _numbering[1] = regions[1].clone();
    _numbering[3] = regions[3].clone();
    _numbering[5] = regions[5].clone();
    _numbering[7] = regions[7].clone();

    // Region 1 (index 0): insertions placed at Chothia position 8.
    let insertions = regions[0].iter().filter(|r| r.ins != " ").count();
    if insertions > 0 {
        let start = regions[0][0].num;
        let length = regions[0].len();
        let mut ann: Vec<(i32, &str)> = Vec::new();
        for p in start..9 {
            ann.push((p, " "));
        }
        for i in 0..insertions {
            ann.push((8, alpha(i as i32)));
        }
        ann.push((9, " "));
        _numbering[0] = apply_annotations(&ann[..length], &regions[0]);
    } else {
        _numbering[0] = regions[0].clone();
    }

    // CDR1 (index 2): insertions onto 31.
    _numbering[2] = cdrh1_chothia_martin(&regions[2]);

    // CDR2 (index 4): insertions onto 52.
    _numbering[4] = cdrh2_50_57(&regions[4]);

    // FW3 (index 5): insertions on 72. The else-branch in Python reassigns
    // _numbering[4] (a dead branch) — so when there are no FW3 insertions, FW3
    // region (_numbering[5]) stays as _regions[5], already set above.
    {
        let length = regions[5].len();
        let insertions = length.saturating_sub(35);
        if insertions > 0 {
            let mut ann: Vec<(i32, &str)> = (58..73).map(|p| (p, " ")).collect();
            for i in 0..insertions {
                ann.push((72, alpha(i as i32)));
            }
            ann.extend((73..93).map(|p| (p, " ")));
            _numbering[5] = apply_annotations(&ann, &regions[5]);
        }
        // else: Python sets _numbering[4]=_regions[4] (no effect on FW3 output).
    }

    // CDR3 (index 6): insertions onto 100.
    if regions[6].len() > 36 {
        return Ok(Numbered { residues: vec![], start: startindex, end: endindex });
    }
    let ann = super::get_cdr3_annotations(regions[6].len(), "chothia", "heavy")?;
    _numbering[6] = apply_annotations(&ann, &regions[6]);

    Ok(Numbered { residues: super::gap_missing(&_numbering), start: startindex, end: endindex })
}

// ---------------------------------------------------------------------------
// Martin — light (delegates to chothia light)
// ---------------------------------------------------------------------------

pub(crate) fn number_martin_light(state_vector: &[State], sequence: &[u8]) -> CResult<Numbered> {
    number_chothia_light(state_vector, sequence)
}

// ---------------------------------------------------------------------------
// Wolfguy — heavy
// ---------------------------------------------------------------------------

pub(crate) fn number_wolfguy_heavy(state_vector: &[State], sequence: &[u8]) -> CResult<Numbered> {
    let state_string = b"XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
    let region_string = "11111111111111111111111111222222222222223333333333333344444444444444444444555555555555555555555555555555666666666666677777777777";
    let dict = [('1', 0), ('2', 1), ('3', 2), ('4', 3), ('5', 4), ('6', 5), ('7', 6)];
    let region_of_state = super::region_map(region_string, &dict);
    let mut rels = [100i32, 124, 160, 196, 226, 244, 283];
    let exclude_deletions = [1usize, 3, 5];

    let (regions, startindex, endindex) = super::number_regions(
        sequence, state_vector, state_string, &region_of_state, &mut rels, 7, &exclude_deletions,
    )?;

    let mut _numbering: Vec<Vec<Res>> = vec![Vec::new(); 7];
    _numbering[0] = regions[0].clone();
    _numbering[2] = regions[2].clone();
    _numbering[4] = regions[4].clone();
    _numbering[6] = regions[6].clone();

    // CDRH1: delete symmetrically about 177, right first.
    {
        let mut ordered: Vec<i32> = vec![151];
        for (p1, p2) in (152..176).zip((176..200).rev()) {
            ordered.push(p1);
            ordered.push(p2);
        }
        _numbering[1] = wolfguy_symmetric(&ordered, &regions[1]);
    }

    // CDRH2: delete symmetrically about 271, right first; then right from 288.
    {
        let mut ordered: Vec<i32> = vec![251];
        for (p1, p2) in (252..271).zip((272..291).rev()) {
            ordered.push(p1);
            ordered.push(p2);
        }
        ordered.push(271);
        let mut prefix: Vec<i32> = (291..300).rev().collect();
        prefix.extend(ordered);
        _numbering[3] = wolfguy_symmetric(&prefix, &regions[3]);
    }

    // CDRH3: delete symmetrically about 374, right first.
    {
        let mut ordered: Vec<i32> = Vec::new();
        for (p1, p2) in (356..374).zip((374..392).rev()) {
            ordered.push(p1);
            ordered.push(p2);
        }
        let mut v: Vec<i32> = vec![354, 394, 355, 393, 392];
        v.extend(ordered);
        let mut full: Vec<i32> = vec![331, 332, 399, 398, 351, 352, 397, 353, 396, 395];
        full.extend(v);
        if regions[5].len() > full.len() {
            return Ok(Numbered { residues: vec![], start: startindex, end: endindex });
        }
        _numbering[5] = wolfguy_symmetric(&full, &regions[5]);
    }

    Ok(Numbered { residues: _numbering.into_iter().flatten().collect(), start: startindex, end: endindex })
}

/// Wolfguy CDR symmetric gapping: take the first `length` ordered deletions,
/// sort them, and annotate (all with a blank insertion code).
fn wolfguy_symmetric(ordered_deletions: &[i32], region: &[Res]) -> Vec<Res> {
    let length = region.len();
    let mut positions: Vec<i32> = ordered_deletions[..length].to_vec();
    positions.sort();
    (0..length).map(|i| Res::new(positions[i], " ", region[i].aa)).collect()
}

// ---------------------------------------------------------------------------
// Wolfguy — light
// ---------------------------------------------------------------------------

pub(crate) fn number_wolfguy_light(state_vector: &[State], sequence: &[u8]) -> CResult<Numbered> {
    let state_string = b"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
    let region_string = "1111111AAABBBBBBBBBBBBB222222222222222223333333333333334444444444444455555555555666677777777777777777777888888888888899999999999";
    let dict = [
        ('1', 0), ('A', 1), ('B', 2), ('2', 3), ('3', 4), ('4', 5), ('5', 6), ('6', 7), ('7', 8),
        ('8', 9), ('9', 10),
    ];
    let region_of_state = super::region_map(region_string, &dict);
    let mut rels = [500i32, 500, 500, 527, 560, 595, 631, 630, 630, 646, 683];
    let exclude_deletions = [1usize, 3, 5, 7, 9];

    let (regions, startindex, endindex) = super::number_regions(
        sequence, state_vector, state_string, &region_of_state, &mut rels, 11, &exclude_deletions,
    )?;

    let mut _numbering: Vec<Vec<Res>> = vec![Vec::new(); 11];
    _numbering[0] = regions[0].clone();
    _numbering[2] = regions[2].clone();
    _numbering[4] = regions[4].clone();
    _numbering[6] = regions[6].clone();
    _numbering[8] = regions[8].clone();
    _numbering[10] = regions[10].clone();

    // Region A (index 1): gaps go 508/509/510, insertions on 508.
    {
        let length = regions[1].len();
        let head = [(510i32, " "), (509, " "), (508, " ")];
        let mut ann: Vec<(i32, &str)> = head[..length.min(3)].to_vec();
        let n_ins = (length as i32 - 3).max(0) as usize;
        for i in 0..n_ins {
            ann.push((508, alpha(i as i32)));
        }
        ann.sort();
        _numbering[1] = apply_annotations(&ann, &regions[1]);
    }

    // CDRL1 (index 3): canonical-class numbering.
    {
        let length = regions[3].len();
        let positions = super::get_wolfguy_l1(&regions[3], length)?;
        _numbering[3] = (0..length)
            .map(|i| Res::new(positions[i], " ", regions[3][i].aa))
            .collect();
    }

    // CDRL2 (index 5): delete about 673, then right from 694; keep 651 last.
    {
        let mut middle: Vec<i32> = Vec::new();
        for (p1, p2) in (652..673).zip((673..695).rev()) {
            middle.push(p2);
            middle.push(p1);
        }
        let mut ordered: Vec<i32> = vec![651];
        ordered.extend((695..700).rev());
        ordered.extend(middle);
        ordered.push(673);
        _numbering[5] = wolfguy_symmetric(&ordered, &regions[5]);
    }

    // Region 6 (index 7): indel placement on 711..714, insertions on 714.
    {
        let length = regions[7].len();
        let insertions = (length as i32 - 4).max(0) as usize;
        let head = [(711i32, " "), (712, " "), (713, " "), (714, " ")];
        let mut ann: Vec<(i32, &str)> = head[..length.min(4)].to_vec();
        for i in 0..insertions {
            ann.push((714, alpha(i as i32)));
        }
        _numbering[7] = apply_annotations(&ann, &regions[7]);
    }

    // CDRL3 (index 9): delete symmetrically about 775, right first; then 798,799.
    {
        let mut ordered: Vec<i32> = Vec::new();
        for (p1, p2) in (751..775).zip((776..800).rev()) {
            ordered.push(p1);
            ordered.push(p2);
        }
        ordered.push(775);
        if regions[9].len() > ordered.len() {
            return Ok(Numbered { residues: vec![], start: startindex, end: endindex });
        }
        _numbering[9] = wolfguy_symmetric(&ordered, &regions[9]);
    }

    Ok(Numbered { residues: _numbering.into_iter().flatten().collect(), start: startindex, end: endindex })
}

// ---------------------------------------------------------------------------
// Aho
// ---------------------------------------------------------------------------

pub(crate) fn number_aho(
    state_vector: &[State],
    sequence: &[u8],
    chain_type: &str,
) -> CResult<Numbered> {
    let state_string = b"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
    let region_string = "BBBBBBBBBBCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDEEEEEEEEEEEEEEEFFFFFFFFFFFFFFFFFFFFHHHHHHHHHHHHHHHHIIIIIIIIIIIIIJJJJJJJJJJJJJKKKKKKKKKKK";
    let dict = [
        ('A', 0), ('B', 1), ('C', 2), ('D', 3), ('E', 4), ('F', 5), ('G', 6), ('H', 7), ('I', 8),
        ('J', 9), ('K', 10),
    ];
    let region_of_state = super::region_map(region_string, &dict);
    let mut rels = [0i32, 0, 0, 0, 2, 2, 2, 2, 2, 2, 21];
    let exclude_deletions = [1usize, 3, 4, 5, 7, 9];

    let (regions, startindex, endindex) = super::number_regions(
        sequence, state_vector, state_string, &region_of_state, &mut rels, 11, &exclude_deletions,
    )?;
    let mut endindex = endindex;

    let mut _numbering: Vec<Vec<Res>> = vec![Vec::new(); 11];
    _numbering[0] = regions[0].clone();
    _numbering[1] = regions[1].clone();
    _numbering[2] = regions[2].clone();
    _numbering[4] = regions[4].clone();
    _numbering[6] = regions[6].clone();
    _numbering[8] = regions[8].clone();
    _numbering[9] = regions[9].clone();
    _numbering[10] = regions[10].clone();

    // Move the indel in FW1 onto 8.
    {
        let length = regions[1].len();
        if length > 0 {
            let start = regions[1][0].num;
            let stretch_len = 10 - (start - 1);
            let ann: Vec<(i32, &str)> = if (length as i32) > stretch_len {
                // Insertions present. Place on 8.
                let mut v: Vec<(i32, &str)> = (start..9).map(|p| (p, " ")).collect();
                let n_ins = length as i32 - stretch_len;
                for i in 0..n_ins {
                    v.push((8, alpha(i)));
                }
                v.push((9, " "));
                v.push((10, " "));
                v
            } else {
                // ordered_deletions = [(8," ")] + [(p," ") for p in range(start,11) if p != 8]
                let mut ordered: Vec<(i32, &str)> = vec![(8, " ")];
                ordered.extend((start..11).filter(|&p| p != 8).map(|p| (p, " ")));
                let drop = (stretch_len - length as i32).max(0) as usize;
                let mut v = ordered[drop..].to_vec();
                v.sort();
                v
            };
            _numbering[1] = apply_annotations(&ann, &regions[1]);
        }
    }

    // CDR1 (index 3): chain-type-dependent gap order; insertions on 36.
    {
        let order = aho_cdr1_order(chain_type);
        let length = regions[3].len();
        let drop = (18i32 - length as i32).max(0) as usize;
        let mut ann: Vec<(i32, &str)> = order[drop..].iter().map(|&p| (p, " ")).collect();
        ann.sort();
        let insertions = (length as i32 - 18).max(0) as usize;
        if insertions > 26 {
            return Ok(Numbered { residues: vec![], start: startindex, end: endindex });
        } else if insertions > 0 {
            let insertat = ann.iter().position(|&x| x == (36, " ")).unwrap() + 1;
            if insertat != 12 {
                return Err(assertion!("AHo numbering failed"));
            }
            let mut new_ann = ann[..insertat].to_vec();
            for a in 0..insertions {
                new_ann.push((36, alpha(a as i32)));
            }
            new_ann.extend_from_slice(&ann[insertat..]);
            ann = new_ann;
        }
        _numbering[3] = apply_annotations(&ann, &regions[3]);
    }

    // CDR2 (index 5): gaps symmetric at 63 (A chain different); insertions on 63.
    {
        let order: Vec<i32> = if chain_type == "A" {
            vec![74, 73, 63, 62, 64, 61, 65, 60, 66, 59, 67, 58, 68, 69, 70, 71, 72, 75, 76, 77]
        } else {
            vec![63, 62, 64, 61, 65, 60, 66, 59, 67, 58, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77]
        };
        let length = regions[5].len();
        let drop = (20i32 - length as i32).max(0) as usize;
        let mut ann: Vec<(i32, &str)> = order[drop..].iter().map(|&p| (p, " ")).collect();
        ann.sort();
        let insertions = (length as i32 - 20).max(0) as usize;
        if insertions > 26 {
            return Ok(Numbered { residues: vec![], start: startindex, end: endindex });
        } else if insertions > 0 {
            let insertat = ann.iter().position(|&x| x == (63, " ")).unwrap() + 1;
            if insertat != 6 {
                return Err(assertion!("AHo numbering failed"));
            }
            let mut new_ann = ann[..insertat].to_vec();
            for a in 0..insertions {
                new_ann.push((63, alpha(a as i32)));
            }
            new_ann.extend_from_slice(&ann[insertat..]);
            ann = new_ann;
        }
        _numbering[5] = apply_annotations(&ann, &regions[5]);
    }

    // FW3 (index 7): deletions onto 86 then 85; insertions on 85.
    {
        let order: [i32; 16] = [86, 85, 87, 84, 88, 83, 89, 82, 90, 81, 91, 80, 92, 79, 93, 78];
        let length = regions[7].len();
        let drop = (16i32 - length as i32).max(0) as usize;
        let mut ann: Vec<(i32, &str)> = order[drop..].iter().map(|&p| (p, " ")).collect();
        ann.sort();
        let insertions = (length as i32 - 16).max(0) as usize;
        if insertions > 26 {
            return Ok(Numbered { residues: vec![], start: startindex, end: endindex });
        } else if insertions > 0 {
            let insertat = ann.iter().position(|&x| x == (85, " ")).unwrap() + 1;
            if insertat != 8 {
                return Err(assertion!("AHo numbering failed"));
            }
            let mut new_ann = ann[..insertat].to_vec();
            for a in 0..insertions {
                new_ann.push((85, alpha(a as i32)));
            }
            new_ann.extend_from_slice(&ann[insertat..]);
            ann = new_ann;
        }
        _numbering[7] = apply_annotations(&ann, &regions[7]);
    }

    // CDR3 (index 9): deletions on 123; insertions on 123.
    {
        let order: [i32; 32] = [
            123, 124, 122, 125, 121, 126, 120, 127, 119, 128, 118, 129, 117, 130, 116, 131, 115,
            132, 114, 133, 113, 134, 112, 135, 111, 136, 110, 137, 109, 138, 108, 107,
        ];
        let length = regions[9].len();
        let drop = (32i32 - length as i32).max(0) as usize;
        let mut ann: Vec<(i32, &str)> = order[drop..].iter().map(|&p| (p, " ")).collect();
        ann.sort();
        let insertions = (length as i32 - 32).max(0) as usize;
        if insertions > 26 {
            return Ok(Numbered { residues: vec![], start: startindex, end: endindex });
        } else if insertions > 0 {
            let insertat = ann.iter().position(|&x| x == (123, " ")).unwrap() + 1;
            if insertat != 17 {
                return Err(assertion!("AHo numbering failed"));
            }
            let mut new_ann = ann[..insertat].to_vec();
            for a in 0..insertions {
                new_ann.push((123, alpha(a as i32)));
            }
            new_ann.extend_from_slice(&ann[insertat..]);
            ann = new_ann;
        }
        _numbering[9] = apply_annotations(&ann, &regions[9]);
    }

    // AHo includes one extra light-chain position (149) than IMGT.
    let mut numbering = super::gap_missing(&_numbering);
    if !numbering.is_empty() {
        let last = *numbering.last().unwrap();
        if last.num == 148
            && last.ins == " "
            && last.aa != b'-'
            && endindex.map(|e| e + 1 < sequence.len()).unwrap_or(false)
        {
            let e = endindex.unwrap();
            numbering.push(Res::new(149, " ", sequence[e + 1]));
            endindex = Some(e + 1);
        }
    }

    Ok(Numbered { residues: numbering, start: startindex, end: endindex })
}

/// AHo CDR1 gap order by chain type.
fn aho_cdr1_order(chain_type: &str) -> [i32; 18] {
    match chain_type {
        "L" => [28, 36, 35, 37, 34, 38, 27, 29, 33, 39, 32, 40, 26, 30, 25, 31, 41, 42],
        "K" => [28, 27, 36, 35, 37, 34, 38, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
        "H" => [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
        "A" => [28, 36, 35, 37, 34, 38, 33, 39, 27, 32, 40, 29, 26, 30, 25, 31, 41, 42],
        "B" => [28, 36, 35, 37, 34, 38, 33, 39, 27, 32, 40, 29, 26, 30, 25, 31, 41, 42],
        "D" => [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
        "G" => [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
        _ => [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
    }
}
