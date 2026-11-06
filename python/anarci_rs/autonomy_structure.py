"""Frozen parent C-alpha proximity evidence, not structural or binding validation.

Native parsing retains observed CA residues, not author residue identifiers or
SEQRES completeness. Exact substring correspondence is mandatory: no gapped
alignment, chain alias guessing, coordinate imputation or candidate remapping.
"""
import hashlib
import json
import math

from . import parse_structure
from .species_adaptation import number_domain, _numbering


_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _sequence(value, name):
    if (not isinstance(value, str) or not value or len(value) > 2000
            or any(aa not in _AA for aa in value)):
        raise ValueError(name + " must contain 1..2000 canonical uppercase residues")


def _domain(sequence, domain_format):
    result = number_domain(sequence, domain_format=domain_format)
    if not isinstance(result, dict):
        raise ValueError("Malformed numbered domain")
    if result.get("domain_count") != 1 or result.get("chain") not in (
            {"H"} if domain_format == "VH" else {"K", "L"}):
        raise ValueError("Require exactly one compatible " + domain_format + " domain")
    start, end = result.get("start"), result.get("end")
    if (type(start) is not int or type(end) is not int
            or not 0 <= start <= end < len(sequence)):
        raise ValueError("Malformed domain endpoints")
    numbering = _numbering(result.get("numbering"))
    indices = result.get("sequence_indices")
    if not isinstance(indices, dict) or set(indices) != set(numbering):
        raise ValueError("Malformed numbered sequence index map")
    if (any(type(i) is not int for i in indices.values())
            or sorted(indices.values()) != list(range(start, end + 1))
            or any(sequence[indices[p]] != aa for p, aa in numbering.items())):
        raise ValueError("Numbered domain must map exactly and contiguously to sequence")
    return result


def _coordinates(chain, indices):
    coordinates = []
    for i in indices:
        ca = chain["residues"][i].get("ca")
        if (not isinstance(ca, (list, tuple)) or len(ca) != 3
                or any(type(v) not in (int, float) or not math.isfinite(v) for v in ca)):
            raise ValueError("Missing or nonfinite CA coordinate in chain " + chain["id"])
        coordinates.append(tuple(ca))
    return coordinates


def validate_interface_evidence(evidence, sequence, expected_domain):
    """Validate the complete declaration and rederive its sets; not authentication."""
    if not isinstance(evidence, dict):
        raise ValueError("Interface evidence must be an object")
    try:
        json.dumps(evidence, allow_nan=False)
        if (evidence["schema_version"] != 1 or evidence["query_format"] != "VH"
                or evidence["structure_assessment"] != "parent_CA_proximity_only"
                or evidence["interface_scope"] != "frozen_parent"
                or evidence["query_sequence_sha256"] != hashlib.sha256(sequence.encode()).hexdigest()
                or evidence["domain"] != expected_domain):
            raise ValueError("Interface evidence has incompatible parent identity or residue map")
        roles = evidence["chain_roles"]
        names = [roles["vh_chain_id"], roles["vl_chain_id"], *roles["antigen_chain_ids"]]
        if any(not isinstance(name, str) or not name for name in names) or len(names) != len(set(names)):
            raise ValueError("Interface evidence needs distinct explicit chain roles")
        if not isinstance(roles["antigen_chain_ids"], list):
            raise ValueError("Antigen roles must be a list")
        provenance = evidence["provenance"]
        if (provenance.get("reviewed") is not True or not isinstance(provenance.get("source_id"), str)
                or not provenance["source_id"].strip() or provenance.get("structure_kind") not in {"experimental", "predicted"}
                or provenance.get("assembly_context") not in {"caller_selected_biological_assembly", "asymmetric_unit", "unknown"}):
            raise ValueError("Missing reviewed structure/assembly provenance")
        metrics = evidence["metrics_definition"]
        if metrics.get("units") != "angstrom" or metrics.get("comparison") != "distance <= cutoff":
            raise ValueError("Unsupported proximity metric convention")
        cutoffs = [metrics["contact_cutoff"], metrics["antigen_cutoff"]]
        if any(type(value) not in (int, float) or not math.isfinite(value) or value <= 0 for value in cutoffs):
            raise ValueError("Invalid proximity cutoff")
        rows = evidence["position_evidence"]
        if not isinstance(rows, list) or len(rows) != len(expected_domain["numbering"]):
            raise ValueError("All parent domain positions require proximity evidence")
        seen, chain_indices, interface, antigen = set(), set(), set(), set()
        for row in rows:
            pos = row["position"]
            index = row["chain_residue_index"]
            if pos in seen or pos not in expected_domain["numbering"] or row["sequence_index"] != expected_domain["sequence_indices"][pos]:
                raise ValueError("Proximity evidence has duplicate/unmapped positions")
            if type(index) is not int or index < 0 or index in chain_indices:
                raise ValueError("Proximity evidence has invalid observed-chain indices")
            seen.add(pos); chain_indices.add(index)
            distance, adistance = row["min_vl_ca_distance"], row["min_antigen_ca_distance"]
            if type(distance) not in (int, float) or not math.isfinite(distance) or distance < 0:
                raise ValueError("Missing or invalid VL proximity")
            if roles["antigen_chain_ids"]:
                if type(adistance) not in (int, float) or not math.isfinite(adistance) or adistance < 0:
                    raise ValueError("Missing antigen proximity")
            elif adistance is not None:
                raise ValueError("Antigen proximity cannot be claimed without an antigen chain")
            if distance <= cutoffs[0]: interface.add(pos)
            if adistance is not None and adistance <= cutoffs[1]: antigen.add(pos)
        for key, derived in (("interface_positions", interface), ("antigen_protected_positions", antigen)):
            stated = evidence[key]
            if not isinstance(stated, list) or len(stated) != len(set(stated)) or set(stated) != derived:
                raise ValueError("Declared " + key + " differs from the full proximity evidence")
        if not interface:
            raise ValueError("No measured interface positions")
    except (KeyError, TypeError) as exc:
        raise ValueError("Incomplete interface evidence schema") from exc


def build_interface_evidence(sequence, structure_text, *, vh_chain_id, vl_chain_id,
                             antigen_chain_ids=(), structure_format="pdb",
                             provenance: dict, contact_cutoff=8.0, antigen_cutoff=10.0):
    """Return JSON-compatible, reviewed parent-only IMGT proximity evidence.

    Cutoffs are inclusive Euclidean CA distances in angstroms. ``reviewed=True``
    is a caller assertion of evidence review, NOT permission to mutate every
    contact: the engineering policy must separately review mutable positions.
    All provenance, including assembly/model claims, remains caller-asserted.
    """
    _sequence(sequence, "sequence")
    if not isinstance(structure_text, str) or not structure_text.strip():
        raise ValueError("structure_text must be nonempty text")
    if structure_format not in {"pdb", "cif", "mmcif"}:
        raise ValueError("structure_format must be pdb, cif or mmcif")
    for value in (contact_cutoff, antigen_cutoff):
        if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
            raise ValueError("CA distance cutoffs must be positive finite numbers")
    if not isinstance(provenance, dict):
        raise ValueError("provenance must be a dictionary")
    if (not isinstance(provenance.get("source_id"), str)
            or not provenance["source_id"].strip()
            or provenance.get("structure_kind") not in {"experimental", "predicted"}
            or provenance.get("assembly_context") not in {
                "caller_selected_biological_assembly", "asymmetric_unit", "unknown"}
            or provenance.get("reviewed") is not True):
        raise ValueError("Require source_id, structure_kind, assembly_context and reviewed=True provenance")
    try:
        provenance_copy = json.loads(json.dumps(provenance, allow_nan=False))
    except (TypeError, ValueError) as exc:
        raise ValueError("provenance must be finite JSON data") from exc
    if not isinstance(antigen_chain_ids, (tuple, list)):
        raise ValueError("antigen_chain_ids must be an explicit list or tuple")
    roles = [vh_chain_id, vl_chain_id] + list(antigen_chain_ids)
    if any(not isinstance(c, str) or not c for c in roles) or len(set(roles)) != len(roles):
        raise ValueError("VH, VL and antigen chain IDs must be explicit and distinct")

    domain = _domain(sequence, "VH")
    parsed = parse_structure(structure_text, format=structure_format)
    if not isinstance(parsed, list):
        raise ValueError("Malformed native structure result")
    chains = {}
    for chain in parsed:
        if not isinstance(chain, dict) or not isinstance(chain.get("id"), str):
            raise ValueError("Malformed native chain")
        if chain["id"] in chains:
            raise ValueError("Duplicate native chain ID: " + chain["id"])
        chains[chain["id"]] = chain
    for chain_id in roles:
        if chain_id not in chains:
            raise ValueError("Requested exact native chain ID is absent: " + chain_id)
        chain = chains[chain_id]
        residues, observed = chain.get("residues"), chain.get("sequence")
        if (not isinstance(observed, str) or not observed
                or not isinstance(residues, list) or len(residues) != len(observed)
                or any(not isinstance(r, dict) or r.get("aa") != aa
                       for r, aa in zip(residues, observed))):
            raise ValueError("Native residues do not match observed chain sequence")

    vh, vl = chains[vh_chain_id], chains[vl_chain_id]
    domain_sequence = sequence[domain["start"]:domain["end"] + 1]
    offset = vh["sequence"].find(domain_sequence)
    if offset < 0 or vh["sequence"].find(domain_sequence, offset + 1) >= 0:
        raise ValueError("VH domain requires a unique exact observed-chain substring; gaps cannot be aligned")
    _sequence(vl["sequence"], "observed VL sequence")
    light_domain = _domain(vl["sequence"], "VL")
    vl_ca = _coordinates(vl, range(light_domain["start"], light_domain["end"] + 1))
    antigen_ca = []
    for chain_id in antigen_chain_ids:
        chain = chains[chain_id]
        antigen_ca.extend(_coordinates(chain, range(len(chain["sequence"]))))
    positions = sorted(domain["sequence_indices"], key=domain["sequence_indices"].get)
    chain_indices = [offset + domain["sequence_indices"][p] - domain["start"] for p in positions]
    vh_ca = _coordinates(vh, chain_indices)
    evidence, interface, antigen = [], [], []
    for position, chain_index, ca in zip(positions, chain_indices, vh_ca):
        light_distance = min(math.dist(ca, other) for other in vl_ca)
        antigen_distance = min((math.dist(ca, other) for other in antigen_ca), default=None)
        if not math.isfinite(light_distance) or (antigen_distance is not None and not math.isfinite(antigen_distance)):
            raise ValueError("Nonfinite calculated CA distance")
        if light_distance <= contact_cutoff:
            interface.append(position)
        if antigen_distance is not None and antigen_distance <= antigen_cutoff:
            antigen.append(position)
        evidence.append({"position": position, "sequence_index": domain["sequence_indices"][position],
                         "chain_residue_index": chain_index, "min_vl_ca_distance": light_distance,
                         "min_antigen_ca_distance": antigen_distance})
    if not interface:
        raise ValueError("No VH/VL variable-domain CA interface contacts at the specified cutoff")
    return {
        "schema_version": 1, "query_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "structure_sha256": hashlib.sha256(structure_text.encode()).hexdigest(),
        "query_format": "VH", "structure_format": structure_format,
        "chain_roles": {"vh_chain_id": vh_chain_id, "vl_chain_id": vl_chain_id,
                        "antigen_chain_ids": list(antigen_chain_ids)},
        "domain": {"start": domain["start"], "end": domain["end"],
                   "numbering": dict(domain["numbering"]), "sequence_indices": dict(domain["sequence_indices"])},
        "interface_positions": interface, "antigen_protected_positions": antigen,
        "position_evidence": evidence, "provenance": provenance_copy,
        "provenance_basis": "caller_asserted_not_independently_verified",
        "structure_assessment": "parent_CA_proximity_only", "interface_scope": "frozen_parent",
        "experimental_validation": False,
        "metrics_definition": {"metric": "minimum Euclidean C-alpha proximity, not SASA",
                               "units": "angstrom", "comparison": "distance <= cutoff",
                               "contact_cutoff": contact_cutoff, "antigen_cutoff": antigen_cutoff,
                               "vl_scope": "numbered VL variable domain only",
                               "antigen_scope": "all observed residues of explicitly selected antigen chains"},
        "coordinate_convention": "zero-based input sequence and native observed chain indices; not author residue IDs",
        "unknowns": (["antigen proximity unknown: no antigen chains supplied"] if not antigen_chain_ids else []) + [
            "Unobserved VL/antigen residues and original author numbering are not available from the native parser.",
            "Biological assembly, side-chain contacts, solvent accessibility and candidate structural effects are not established.",
            "Autonomy, solubility, stability, binding and specificity require independent experiments."],
        "parser_limits": [
            "Native CA-only observed residues; no missing-residue reconstruction or coordinate imputation.",
            "PDB stops at first ENDMDL, accepts blank/A alternate locations, deduplicates residue IDs.",
            "mmCIF uses first encountered CA model, accepts ./?/A alternate locations, deduplicates residue IDs.",
            "Native parser normalizes MSE/SEC/PYL to M/C/K and omits unrecognized residue names."]}
