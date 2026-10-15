//! Shared structural types: a chain as its ordered residues, each with a Cα coordinate.
//! Used by the PDB/mmCIF readers (coordinate path) and the pairing inference.

use crate::seqio::Record;

#[derive(Clone, Copy)]
pub struct ResCoord {
    pub aa: u8,
    /// Cα coordinate; `[NaN; 3]` if the file lacked/garbled it (skipped in geometry).
    pub ca: [f32; 3],
}

pub struct StructChain {
    pub id: String,
    pub residues: Vec<ResCoord>,
}

impl StructChain {
    pub fn seq(&self) -> Vec<u8> {
        self.residues.iter().map(|r| r.aa).collect()
    }
    pub fn to_record(&self) -> Record {
        Record { id: self.id.clone(), desc: String::new(), seq: self.seq() }
    }
}

/// Squared Euclidean distance (avoids a sqrt in the O(n·m) contact scan).
pub fn dist2(a: [f32; 3], b: [f32; 3]) -> f32 {
    let (dx, dy, dz) = (a[0] - b[0], a[1] - b[1], a[2] - b[2]);
    dx * dx + dy * dy + dz * dz
}

/// Parse a fixed-width PDB coordinate field (trimmed ASCII float); NaN on failure.
pub fn parse_coord(b: &[u8]) -> f32 {
    std::str::from_utf8(b).ok().and_then(|s| s.trim().parse::<f32>().ok()).unwrap_or(f32::NAN)
}
