//! Infer VH/VL pairing and antibody format from 3D structure.
//!
//! Pairing is a spatial question: two V domains form an Fv when their framework
//! interfaces are in physical contact. We count inter-domain Cα contacts (< 8 Å),
//! then greedily match each heavy domain to the light domain it contacts most (a
//! real Fv interface has many contacts; incidental proximity has few). The topology
//! of the resulting graph — intra- vs inter-chain, cross-paired, unpaired — names the
//! format: scFv, tandem scFv, diabody, VHH, and mixtures thereof.

use crate::structure::{dist2, StructChain};
use serde::Serialize;
use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

const ANTIGEN_A: f32 = 10.0; // Cα-Cα paratope↔antigen distance (Å; side chains bridge the gap)
const MIN_ANTIGEN_CONTACTS: u32 = 5;

/// Chains, Fvs, and unpaired domains gathered into one assembly (grouping accumulator).
type Group = (Vec<String>, Vec<Fv>, Vec<UnpairedDomain>);

const CONTACT_A: f32 = 8.0; // Cα-Cα contact distance (Å)
const MIN_CONTACTS: u32 = 8; // an Fv interface clears this; incidental contact doesn't

/// One numbered V domain with the Cα coordinates of its residues.
pub struct VDomain {
    pub id: String,
    pub chain: String,
    pub chain_type: String, // H / K / L (or TCR B/A/D/G)
    pub species: String,
    pub v_gene: Option<String>,
    pub j_gene: Option<String>,
    pub cas: Vec<[f32; 3]>,
    /// Cα of the CDR residues (IMGT 27–38 / 56–65 / 105–117) — the paratope, used to
    /// find which chain(s) the domain actually binds.
    pub paratope: Vec<[f32; 3]>,
}

/// A chain the Fv/domain binds via its paratope. Usually the antigen; `is_antibody`
/// marks the anti-idiotype case (the binding partner is itself an antibody V domain).
#[derive(Serialize)]
pub struct Antigen {
    pub chain: String,
    pub contacts: u32,
    pub is_antibody: bool,
}

impl VDomain {
    fn is_heavy(&self) -> bool {
        matches!(self.chain_type.as_str(), "H" | "B" | "D")
    }
}

/// A β–δ chain-type pair. The dual TRAV/DV genes make some αβ-TCR α chains number as δ
/// (chain_type "D"), landing them on the heavy side beside their β partner — so a real αβ Fv
/// reads as two heavies and never pairs. β only ever pairs α, so a β in Cα contact with a "δ" is
/// an αβ TCR whose α was mistyped: treat it as an Fv interface. (Real γδ pairs G–D, already
/// heavy–light, and is untouched. Empirically δ→α mistyping never occurs, so no other exception.)
fn beta_delta(a: &VDomain, b: &VDomain) -> bool {
    matches!(
        (a.chain_type.as_str(), b.chain_type.as_str()),
        ("B", "D") | ("D", "B")
    )
}

#[derive(Serialize)]
pub struct Fv {
    pub heavy: String,
    pub light: String,
    pub species: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub v_gene: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub j_gene: Option<String>,
    pub contacts: u32,
    pub intra_chain: bool,
    /// Chain(s) this Fv binds (usually one antigen; ≥2 for epitopes spanning subunits).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub antigen: Vec<Antigen>,
}

/// An independent structural unit: the chains connected (directly or via Fvs) into one
/// assembly, with its Fvs, unpaired domains, and classified format. Separates
/// crystallographic copies and genuinely mixed complexes (e.g. diabody + VHH).
#[derive(Serialize)]
pub struct Assembly {
    pub format: String,
    pub chains: Vec<String>,
    pub fvs: Vec<Fv>,
    pub unpaired: Vec<UnpairedDomain>,
}

#[derive(Serialize)]
pub struct Pairing {
    /// Aggregate over assemblies, e.g. "2× diabody" or "diabody + VHH".
    pub format: String,
    pub assemblies: Vec<Assembly>,
    pub n_domains: usize,
}

#[derive(Serialize)]
pub struct UnpairedDomain {
    pub id: String,
    pub chain: String,
    pub species: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub v_gene: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub j_gene: Option<String>,
    pub kind: String, // "VHH (unpaired VH)" | "orphan VL"
    /// Chain(s) this domain binds (a VHH/nanobody binds an antigen directly).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub antigen: Vec<Antigen>,
}

fn contacts(a: &VDomain, b: &VDomain) -> u32 {
    let t2 = CONTACT_A * CONTACT_A;
    let mut n = 0;
    for &p in &a.cas {
        if p[0].is_nan() {
            continue;
        }
        for &q in &b.cas {
            if !q[0].is_nan() && dist2(p, q) < t2 {
                n += 1;
            }
        }
    }
    n
}

/// Infer Fv pairing and format from the numbered domains.
pub fn infer(domains: &[VDomain]) -> Pairing {
    // candidate edges: heavy–light pairs (a VH-VL Fv is one heavy + one light), plus the β–δ
    // TCR exception, scored by contact count. Sorted strongest-first for greedy matching.
    let mut edges: Vec<(usize, usize, u32)> = Vec::new();
    for i in 0..domains.len() {
        for j in (i + 1)..domains.len() {
            if domains[i].is_heavy() == domains[j].is_heavy()
                && !beta_delta(&domains[i], &domains[j])
            {
                continue; // same class (both heavy / both light) and not a β–δ TCR pair
            }
            let c = contacts(&domains[i], &domains[j]);
            if c >= MIN_CONTACTS {
                edges.push((i, j, c));
            }
        }
    }
    edges.sort_by_key(|e| Reverse(e.2)); // strongest interface first

    let mut used = vec![false; domains.len()];
    let mut fvs = Vec::new();
    for (i, j, c) in edges {
        if used[i] || used[j] {
            continue;
        }
        used[i] = true;
        used[j] = true;
        // heavy analog: the heavy-class domain, or β when both read heavy (β–δ TCR pair).
        let i_heavy = if domains[i].is_heavy() != domains[j].is_heavy() {
            domains[i].is_heavy()
        } else {
            domains[i].chain_type == "B"
        };
        let (h, l) = if i_heavy { (i, j) } else { (j, i) };
        fvs.push(Fv {
            heavy: domains[h].id.clone(),
            light: domains[l].id.clone(),
            species: domains[h].species.clone(),
            v_gene: domains[h].v_gene.clone(),
            j_gene: domains[h].j_gene.clone(),
            contacts: c,
            intra_chain: domains[h].chain == domains[l].chain,
            antigen: Vec::new(),
        });
    }

    let unpaired: Vec<UnpairedDomain> = domains
        .iter()
        .enumerate()
        .filter(|(k, _)| !used[*k])
        .map(|(_, d)| UnpairedDomain {
            id: d.id.clone(),
            chain: d.chain.clone(),
            species: d.species.clone(),
            v_gene: d.v_gene.clone(),
            j_gene: d.j_gene.clone(),
            kind: if d.is_heavy() { "VHH (unpaired VH)".into() } else { "orphan VL".into() },
            antigen: Vec::new(),
        })
        .collect();

    // group chains into independent assemblies: union chains linked by an inter-chain Fv
    let chain_of: HashMap<&str, &str> = domains.iter().map(|d| (d.id.as_str(), d.chain.as_str())).collect();
    let chains: Vec<String> = domains.iter().map(|d| d.chain.clone()).collect::<BTreeSet<_>>().into_iter().collect();
    let cidx: HashMap<&str, usize> = chains.iter().enumerate().map(|(i, c)| (c.as_str(), i)).collect();
    let mut uf: Vec<usize> = (0..chains.len()).collect();
    for f in &fvs {
        if !f.intra_chain {
            let (a, b) = (cidx[chain_of[f.heavy.as_str()]], cidx[chain_of[f.light.as_str()]]);
            let (ra, rb) = (uf_find(&mut uf, a), uf_find(&mut uf, b));
            uf[ra] = rb;
        }
    }
    let mut groups: BTreeMap<usize, Group> = BTreeMap::new();
    for (i, c) in chains.iter().enumerate() {
        let r = uf_find(&mut uf, i);
        groups.entry(r).or_default().0.push(c.clone());
    }
    for f in fvs {
        let r = uf_find(&mut uf, cidx[chain_of[f.heavy.as_str()]]);
        groups.get_mut(&r).unwrap().1.push(f);
    }
    for u in unpaired {
        let r = uf_find(&mut uf, cidx[u.chain.as_str()]);
        groups.get_mut(&r).unwrap().2.push(u);
    }
    let assemblies: Vec<Assembly> = groups
        .into_values()
        .map(|(chs, fvs, unp)| Assembly {
            format: classify_assembly(&fvs, &unp, chs.len(), &chain_of),
            chains: chs,
            fvs,
            unpaired: unp,
        })
        .collect();

    Pairing { format: aggregate(&assemblies), assemblies, n_domains: domains.len() }
}

fn uf_find(uf: &mut [usize], x: usize) -> usize {
    if uf[x] != x {
        let r = uf_find(uf, uf[x]);
        uf[x] = r;
    }
    uf[x]
}

/// Annotate each Fv (and unpaired domain) with the chain(s) its paratope binds. Candidate
/// antigens are chains carrying no V domain; a paratope contacting two of them flags a
/// multi-subunit epitope (e.g. integrin RGD across α/β, or complex-specific antibodies).
pub fn assign_antigens(p: &mut Pairing, domains: &[VDomain], chains: &[StructChain]) {
    let ab_chains: HashSet<&str> = domains.iter().map(|d| d.chain.as_str()).collect();
    let paratope_of = |id: &str| -> &[[f32; 3]] {
        domains.iter().find(|d| d.id == id).map(|d| d.paratope.as_slice()).unwrap_or(&[])
    };
    let chain_of = |id: &str| -> String {
        domains.iter().find(|d| d.id == id).map(|d| d.chain.clone()).unwrap_or_default()
    };
    // chains whose residues the paratope contacts, excluding the domain's OWN chains (its Fv
    // partner is not its antigen). A contacted *antibody* chain is the anti-idiotype case.
    let find = |paratope: &[[f32; 3]], own: &HashSet<String>| -> Vec<Antigen> {
        let t2 = ANTIGEN_A * ANTIGEN_A;
        let mut ags: Vec<Antigen> = chains
            .iter()
            .filter(|c| !own.contains(c.id.as_str()))
            .filter_map(|c| {
                let mut n = 0u32;
                for &pt in paratope {
                    if pt[0].is_nan() {
                        continue;
                    }
                    for r in &c.residues {
                        if !r.ca[0].is_nan() && dist2(pt, r.ca) < t2 {
                            n += 1;
                        }
                    }
                }
                (n >= MIN_ANTIGEN_CONTACTS).then(|| Antigen {
                    chain: c.id.clone(),
                    contacts: n,
                    is_antibody: ab_chains.contains(c.id.as_str()),
                })
            })
            .collect();
        ags.sort_by_key(|a| Reverse(a.contacts));
        ags
    };
    for a in &mut p.assemblies {
        for f in &mut a.fvs {
            let own: HashSet<String> = [chain_of(&f.heavy), chain_of(&f.light)].into_iter().collect();
            let mut para = paratope_of(&f.heavy).to_vec();
            para.extend_from_slice(paratope_of(&f.light));
            f.antigen = find(&para, &own);
        }
        for u in &mut a.unpaired {
            let own: HashSet<String> = [u.chain.clone()].into_iter().collect();
            u.antigen = find(paratope_of(&u.id), &own);
        }
    }
}

/// Classify one assembly's sub-topology, folding in any unpaired VHH / orphan domains.
fn classify_assembly(fvs: &[Fv], unpaired: &[UnpairedDomain], nchains: usize, chain_of: &HashMap<&str, &str>) -> String {
    let intra = fvs.iter().filter(|f| f.intra_chain).count();
    let inter = fvs.len() - intra;
    let vhh = unpaired.iter().filter(|u| u.kind.starts_with("VHH")).count();
    let orphan = unpaired.len() - vhh;
    let base = match (fvs.len(), intra, inter) {
        (0, _, _) => String::new(),
        (1, 1, 0) => "scFv".into(),
        (1, 0, 1) => "Fab/Fv".into(),
        (2, 2, 0) if nchains == 1 => "tandem scFv".into(),
        (2, 0, 2) if is_diabody(fvs, chain_of) => "diabody".into(),
        (2, 0, 2) if is_dual_variable(fvs, chain_of) => dual_variable_name(fvs),
        (3, 0, 3) if is_triabody(fvs, chain_of) => "triabody".into(),
        (n, i, 0) if i == n => format!("tandem ({n}× scFv)"),
        (n, 0, _) => format!("{n}× inter-chain Fv"),
        (n, _, _) => format!("{n} Fv ({intra} intra / {inter} inter)"),
    };
    let mut parts = Vec::new();
    if !base.is_empty() {
        parts.push(base);
    }
    match vhh {
        0 => {}
        1 => parts.push("VHH".into()),
        n => parts.push(format!("{n}× VHH")),
    }
    if orphan > 0 {
        parts.push(format!("{orphan}× orphan VL"));
    }
    if parts.is_empty() {
        "no V domains".into()
    } else {
        parts.join(" + ")
    }
}

/// Two inter-chain Fvs over exactly two chains, cross-paired (the diabody signature:
/// each chain's heavy pairs the OTHER chain's light).
fn is_diabody(fvs: &[Fv], chain_of: &HashMap<&str, &str>) -> bool {
    if fvs.len() != 2 {
        return false;
    }
    let h = |f: &Fv| chain_of[f.heavy.as_str()];
    let l = |f: &Fv| chain_of[f.light.as_str()];
    let chs: HashSet<&str> = [h(&fvs[0]), l(&fvs[0]), h(&fvs[1]), l(&fvs[1])].into_iter().collect();
    chs.len() == 2 && h(&fvs[0]) != l(&fvs[0]) && h(&fvs[1]) != l(&fvs[1]) && h(&fvs[0]) == l(&fvs[1]) && h(&fvs[1]) == l(&fvs[0])
}

fn chain_pairs(fvs: &[Fv], chain_of: &HashMap<&str, &str>) -> Vec<(String, String)> {
    fvs.iter()
        .map(|f| (chain_of[f.heavy.as_str()].to_string(), chain_of[f.light.as_str()].to_string()))
        .collect()
}

/// The domain's index within its chain (the `N` in `id:N`).
fn dom_idx(id: &str) -> usize {
    id.rsplit(':').next().and_then(|s| s.parse().ok()).unwrap_or(0)
}

/// Two inter-chain Fvs where one chain is *always* heavy and the other *always* light
/// (not cross-paired) — the dual-variable topology (DVD-Ig / CODV), two stacked Fvs.
fn is_dual_variable(fvs: &[Fv], chain_of: &HashMap<&str, &str>) -> bool {
    if fvs.len() != 2 {
        return false;
    }
    let c = chain_pairs(fvs, chain_of);
    let chs: HashSet<&str> = c.iter().flat_map(|(h, l)| [h.as_str(), l.as_str()]).collect();
    chs.len() == 2 && c[0].0 == c[1].0 && c[0].1 == c[1].1
}

/// Matched domain order (Hi–Li: outer–outer, inner–inner) = DVD-Ig; crossed = CODV.
fn dual_variable_name(fvs: &[Fv]) -> String {
    if fvs.iter().all(|f| dom_idx(&f.heavy) == dom_idx(&f.light)) {
        "DVD-Ig (dual-variable)".into()
    } else {
        "CODV (crossover dual-variable)".into()
    }
}

/// Three inter-chain Fvs over three chains, each chain heavy in one Fv and light in
/// another (cyclic cross-pairing) — a triabody.
fn is_triabody(fvs: &[Fv], chain_of: &HashMap<&str, &str>) -> bool {
    if fvs.len() != 3 {
        return false;
    }
    let c = chain_pairs(fvs, chain_of);
    let heavies: HashSet<&str> = c.iter().map(|(h, _)| h.as_str()).collect();
    let lights: HashSet<&str> = c.iter().map(|(_, l)| l.as_str()).collect();
    heavies.len() == 3 && lights.len() == 3 && c.iter().all(|(h, l)| h != l)
}

/// Aggregate the per-assembly formats: "diabody" ×2 → "2× diabody"; distinct → joined.
fn aggregate(assemblies: &[Assembly]) -> String {
    if assemblies.len() == 1 {
        return assemblies[0].format.clone();
    }
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for a in assemblies {
        *counts.entry(a.format.as_str()).or_default() += 1;
    }
    let mut items: Vec<(&str, usize)> = counts.into_iter().collect();
    items.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(b.0)));
    items
        .iter()
        .map(|(f, n)| if *n > 1 { format!("{n}× {f}") } else { f.to_string() })
        .collect::<Vec<_>>()
        .join(" + ")
}

/// Human-readable report, grouped by assembly.
pub fn report(p: &Pairing) -> String {
    let plural = if p.assemblies.len() == 1 { "assembly" } else { "assemblies" };
    let mut s = format!("format: {}\n{} V domain(s) in {} {}\n", p.format, p.n_domains, p.assemblies.len(), plural);
    for (i, a) in p.assemblies.iter().enumerate() {
        s += &format!("─ assembly {} [chains {}]: {}\n", i + 1, a.chains.join(","), a.format);
        for f in &a.fvs {
            s += &format!(
                "    Fv [{}]: {} — {}  [{}, {} contacts]\n",
                f.species, f.heavy, f.light,
                if f.intra_chain { "intra-chain" } else { "inter-chain" },
                f.contacts,
            );
            s += &fmt_antigen(&f.antigen);
        }
        for u in &a.unpaired {
            s += &format!("    unpaired [{}]: {} — {}\n", u.species, u.id, u.kind);
            s += &fmt_antigen(&u.antigen);
        }
    }
    s
}

fn fmt_antigen(ags: &[Antigen]) -> String {
    let show = |v: &[&Antigen]| v.iter().map(|a| format!("{} ({} contacts)", a.chain, a.contacts)).collect::<Vec<_>>().join(", ");
    let antigens: Vec<&Antigen> = ags.iter().filter(|a| !a.is_antibody).collect();
    let idiotype: Vec<&Antigen> = ags.iter().filter(|a| a.is_antibody).collect();
    let mut s = String::new();
    if !antigens.is_empty() {
        let tag = if antigens.len() > 1 { "binds (multi-subunit epitope)" } else { "binds antigen" };
        s += &format!("        {}: {}\n", tag, show(&antigens));
    }
    if !idiotype.is_empty() {
        s += &format!("        anti-idiotype (binds antibody): {}\n", show(&idiotype));
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dom(id: &str, chain: &str, ct: &str, at: [f32; 3]) -> VDomain {
        // a small blob of 10 CAs around `at`; paratope = the same blob (for antigen tests)
        let cas: Vec<[f32; 3]> = (0..10).map(|k| [at[0] + k as f32 * 0.5, at[1], at[2]]).collect();
        VDomain { id: id.into(), chain: chain.into(), chain_type: ct.into(), species: "human".into(), v_gene: None, j_gene: None, paratope: cas.clone(), cas }
    }

    #[test]
    fn diabody_cross_pairing() {
        // chain A: VH_A near origin, VL_A far; chain B: VH_B near VL_A, VL_B near VH_A.
        // so VH_A(A) pairs VL_B(B) and VH_B(B) pairs VL_A(A) — the diabody cross.
        let ds = vec![
            dom("A:0", "A", "H", [0.0, 0.0, 0.0]),   // VH_A
            dom("A:1", "A", "K", [0.0, 0.0, 30.0]),  // VL_A (far from VH_A)
            dom("B:0", "B", "H", [0.0, 0.0, 33.0]),  // VH_B near VL_A
            dom("B:1", "B", "K", [0.0, 2.0, 0.0]),   // VL_B near VH_A
        ];
        let p = infer(&ds);
        // one assembly (chains A+B linked), classified diabody, 2 inter-chain Fvs
        assert_eq!(p.format, "diabody");
        assert_eq!(p.assemblies.len(), 1);
        assert_eq!(p.assemblies[0].fvs.len(), 2);
        assert!(p.assemblies[0].fvs.iter().all(|f| !f.intra_chain));
    }

    #[test]
    fn lone_vhh_is_unpaired() {
        let ds = vec![dom("N:0", "N", "H", [0.0, 0.0, 0.0])];
        let p = infer(&ds);
        assert_eq!(p.format, "VHH");
        assert_eq!(p.assemblies.len(), 1);
        assert_eq!(p.assemblies[0].unpaired.len(), 1);
        assert!(p.assemblies[0].fvs.is_empty());
    }

    #[test]
    fn scfv_intra_chain() {
        // one chain, VH and VL in contact -> intra-chain Fv
        let ds = vec![
            dom("S:0", "S", "H", [0.0, 0.0, 0.0]),
            dom("S:1", "S", "L", [0.0, 2.0, 0.0]),
        ];
        let p = infer(&ds);
        assert_eq!(p.format, "scFv");
        assert_eq!(p.assemblies[0].fvs.len(), 1);
        assert!(p.assemblies[0].fvs[0].intra_chain);
    }

    #[test]
    fn tcr_beta_delta_pairs() {
        // αβ TCR whose α numbered as δ (dual TRAV/DV gene): β("B") + "δ"("D") on separate chains
        // must pair into one inter-chain Fv, with β as the heavy analog — not two lone VHHs.
        let ds = vec![
            dom("D:0", "D", "D", [0.0, 0.0, 0.0]), // α mistyped as δ
            dom("E:0", "E", "B", [0.0, 2.0, 0.0]), // β partner, in contact
        ];
        let p = infer(&ds);
        assert_eq!(p.assemblies.len(), 1);
        assert_eq!(p.assemblies[0].fvs.len(), 1);
        assert!(p.assemblies[0].unpaired.is_empty());
        assert_eq!(p.assemblies[0].fvs[0].heavy, "E:0"); // β is heavy
        assert_eq!(p.assemblies[0].fvs[0].light, "D:0");
        // two β with no α/δ partner must NOT pair (β–β is not an Fv)
        let bb = vec![dom("X:0", "X", "B", [0.0, 0.0, 0.0]), dom("Y:0", "Y", "B", [0.0, 2.0, 0.0])];
        assert!(infer(&bb).assemblies.iter().all(|a| a.fvs.is_empty()));
    }

    #[test]
    fn antigen_is_the_contacted_nonantibody_chain() {
        use crate::structure::{ResCoord, StructChain};
        let ds = vec![dom("F:0", "F", "H", [0.0, 0.0, 0.0]), dom("F:1", "F", "L", [0.0, 2.0, 0.0])];
        let ag = |id: &str, at: [f32; 3]| StructChain {
            id: id.into(),
            residues: (0..8).map(|k| ResCoord { aa: b'A', ca: [at[0] + k as f32 * 0.5, at[1], at[2]] }).collect(),
        };
        // chain X sits on the paratope; chain Z is far away (crystal-distant, not antigen)
        let chains = vec![ag("X", [3.0, 1.0, 0.0]), ag("Z", [200.0, 200.0, 200.0])];
        let mut p = infer(&ds);
        assign_antigens(&mut p, &ds, &chains);
        let ag = &p.assemblies[0].fvs[0].antigen;
        assert_eq!(ag.len(), 1);
        assert_eq!(ag[0].chain, "X");
        assert!(!ag[0].is_antibody); // a real (non-antibody) antigen
    }

    #[test]
    fn anti_idiotype_partner_flagged_as_antibody() {
        use crate::structure::{ResCoord, StructChain};
        // two separate scFv antibodies whose paratopes touch each other (idiotype ↔ anti-idiotype)
        let ds = vec![
            dom("P:0", "P", "H", [0.0, 0.0, 0.0]),
            dom("P:1", "P", "L", [0.0, 2.0, 0.0]),
            dom("Q:0", "Q", "H", [0.0, 5.0, 0.0]),
            dom("Q:1", "Q", "L", [0.0, 7.0, 0.0]),
        ];
        // the antibody chains themselves are present as candidate contacts
        let chain = |id: &str, ds: &[&VDomain]| StructChain {
            id: id.into(),
            residues: ds.iter().flat_map(|d| d.cas.iter().map(|&ca| ResCoord { aa: b'A', ca })).collect(),
        };
        let chains = vec![chain("P", &[&ds[0], &ds[1]]), chain("Q", &[&ds[2], &ds[3]])];
        let mut p = infer(&ds);
        assign_antigens(&mut p, &ds, &chains);
        let ags: Vec<&Antigen> = p.assemblies.iter().flat_map(|a| a.fvs.iter()).flat_map(|f| f.antigen.iter()).collect();
        assert!(!ags.is_empty() && ags.iter().all(|a| a.is_antibody), "partner should be flagged anti-idiotype");
    }
}
