"""Exact, explicitly scoped antibody engineering hypotheses, not therapeutic validation."""
import bisect
import copy
import hashlib
import json
import math
import re
import sys
from pathlib import Path

from .minimum_edit import solve_minimum
from .species_reference import validate_profile, score_numbered
from .species_adaptation import AA, CDRS, _position, _protected, _raw_score, _liabilities, number_domain

# Kabat Vernier position sets tabulated by AbNumber from Foote & Winter (1992).
# Reference: common.py blob d494cb3593745ee824727537176da5eec67e1aaf.
# Positions are mapped by actual raw residue index, never copied as IMGT numbers.
VERNIER = {
    "H": frozenset([2, 27, 28, 29, 30, 47, 48, 49, 67, 69, 71, 73, 78, 93, 94, 103]),
    "K": frozenset([2, 4, 35, 36, 46, 47, 48, 49, 64, 66, 68, 69, 71, 98]),
    "L": frozenset([2, 4, 35, 36, 46, 47, 48, 49, 64, 66, 68, 69, 71, 98]),
}
HYDROPHOBIC = frozenset("AVILMFWY")


def _hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _finite(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(name + " must be finite numeric data")
    return float(value)


def _threshold(value, model):
    if not isinstance(value, dict) or set(value) != {"kind", "minimum", "rationale"}:
        raise ValueError("threshold requires exactly kind, minimum and rationale")
    kind = value["kind"]
    if kind not in {"raw_score", "percentile"}:
        raise ValueError("threshold kind must be raw_score or percentile")
    minimum = _finite(value["minimum"], "threshold minimum")
    if not isinstance(value["rationale"], str) or not value["rationale"].strip():
        raise ValueError("State why this research threshold was selected")
    if kind == "percentile" and (not 0 <= minimum <= 100 or not model["calibration_scores"]):
        raise ValueError("Percentile targets need calibration and a minimum in 0..100")
    return {**value, "minimum": minimum, "minimum_binary64": minimum.hex()}


def _objective(model, numbering, threshold):
    raw = _raw_score(model, numbering)
    value = raw
    if threshold["kind"] == "percentile":
        scores = model["calibration_scores"]
        value = 100.0 * (bisect.bisect_left(scores, raw) + 0.5 * (bisect.bisect_right(scores, raw) - bisect.bisect_left(scores, raw))) / len(scores)
    return value >= threshold["minimum"], {"kind": threshold["kind"], "value": value,
                                           "minimum": threshold["minimum"], "raw_score": raw}


def _positions(values, domain, label):
    if not isinstance(values, (list, tuple)):
        raise ValueError(label + " must be a list of IMGT positions")
    for pos in values:
        _position(pos)
        if pos not in domain["numbering"]:
            raise ValueError(label + " position absent from domain: " + pos)
    return set(values)


def _native_map(sequence, scheme="imgt"):
    """Return successful native numbering as data; operational errors propagate."""
    from . import anarci
    numbered, details, _ = anarci([("design", sequence)], scheme=scheme, database="ALL", allowed_species=None, ncpu=1)
    if not isinstance(numbered, list) or len(numbered) != 1 or not isinstance(details, list) or len(details) != 1:
        raise ValueError("Malformed native batch result")
    if numbered[0] is None:
        if details[0] is not None:
            raise ValueError("Malformed native no-hit details")
        return None
    if (not isinstance(numbered[0], list) or not numbered[0] or not isinstance(details[0], list)
            or len(numbered[0]) != len(details[0])):
        raise ValueError("Malformed native domain collections")
    parsed = []
    for record, detail in zip(numbered[0], details[0]):
        if (not isinstance(record, (list, tuple)) or len(record) != 3 or not isinstance(detail, dict)
                or detail.get("chain_type") not in {"H", "K", "L", "A", "B", "G", "D"}):
            raise ValueError("Malformed native domain or chain class")
        positions, start, end = record
        if (not isinstance(positions, (list, tuple)) or type(start) is not int or type(end) is not int
                or not 0 <= start <= end < len(sequence)):
            raise ValueError("Malformed native endpoints")
        mapping, indices, seen_positions = {}, {}, set()
        index = start
        for entry in positions:
            if not isinstance(entry, (list, tuple)) or len(entry) != 2:
                raise ValueError("Malformed native position entry")
            location, residue = entry
            if not isinstance(location, (list, tuple)) or len(location) != 2:
                raise ValueError("Malformed native position label")
            number, insertion = location
            if (type(number) is not int or not isinstance(insertion, str)
                    or re.fullmatch(r"[A-Z]*", insertion.strip()) is None
                    or not isinstance(residue, str) or residue not in AA | {"-"}):
                raise ValueError("Malformed native label/residue types")
            pos = str(number) + insertion.strip()
            if pos in seen_positions:
                raise ValueError("Duplicate native position label")
            seen_positions.add(pos)
            if residue == "-":
                continue
            if index >= len(sequence) or sequence[index] != residue:
                raise ValueError("Malformed native residue correspondence")
            mapping[pos], indices[pos] = residue, index
            index += 1
        if index != end + 1 or not mapping:
            raise ValueError("Malformed native residue count")
        parsed.append({"numbering": mapping, "sequence_indices": indices, "chain": detail["chain_type"],
                       "start": start, "end": end, "domain_count": len(numbered[0])})
    return parsed[0] if len(parsed) == 1 else {"domain_count": len(parsed)}


def _vernier_positions(domain):
    kabat = _native_map(domain["sequence"], "kabat")
    if not kabat or kabat.get("domain_count") != 1 or kabat.get("chain") != domain["chain"]:
        raise ValueError("Cannot establish Kabat-to-IMGT Vernier correspondence")
    if set(kabat["sequence_indices"].values()) != set(domain["sequence_indices"].values()):
        raise ValueError("Kabat and IMGT domain residue coverage differ; Vernier protection is unresolved")
    protected_indices = {index for pos, index in kabat["sequence_indices"].items()
                         if int(re.match(r"\d+", pos)[0]) in VERNIER[domain["chain"]]}
    return {pos for pos, index in domain["sequence_indices"].items() if index in protected_indices}


def template_alphabet(domain, template):
    """Compile ONE explicit human V/J pair; no nearest-gene guess or mosaic union."""
    from .germlines import all_germlines
    if not isinstance(template, dict) or set(template) != {"v_gene", "j_gene"}:
        raise ValueError("template requires exactly one v_gene and one j_gene")
    try:
        v = all_germlines["V"][domain["chain"]]["human"][template["v_gene"]]
        j = all_germlines["J"][domain["chain"]]["human"][template["j_gene"]]
    except KeyError as exc:
        raise ValueError("Unknown human germline template for this chain") from exc
    if len(v) != 128 or len(j) != 128:
        raise ValueError("Expected IMGT-aligned 128-column germline reference")
    alphabet = {}
    for pos, old in domain["numbering"].items():
        base, insertion = _position(pos)
        if insertion or any(low <= base <= high for low, high in CDRS):
            continue
        # FR4 comes from J. Never fill a missing germline coordinate by guessing.
        aa = (j if base >= 118 else v)[base - 1]
        if aa in AA and aa != old:
            alphabet[pos] = [aa]
    return alphabet, {"source": "bundled human IMGT-aligned germlines", **template,
                      "v_sha256": hashlib.sha256(v.encode()).hexdigest(), "j_sha256": hashlib.sha256(j.encode()).hexdigest(),
                      "semantics": "single_V_J_template_substitutions", "likelihood_family_conditioned": False}


def _compile_domain(spec, *, autonomy=False):
    if not isinstance(spec, dict):
        raise ValueError("Each domain specification must be an object")
    supported = {"id", "sequence", "domain_format", "profile", "threshold", "allowed_substitutions", "template",
                 "mutable_positions", "protected_positions", "core_positions", "antigen_positions", "interface_positions",
                 "protection_overrides", "guards", "parent_provenance"}
    if set(spec) - supported:
        raise ValueError("Unknown domain options: " + str(sorted(set(spec) - supported)))
    identifier = spec.get("id")
    if not isinstance(identifier, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]{0,31}", identifier):
        raise ValueError("Domain id must be a short unique identifier without ':'")
    sequence = spec.get("sequence")
    if not isinstance(sequence, str) or len(sequence) > 2000:
        raise ValueError("Supply a canonical sequence of at most2000 residues per domain input")
    domain = number_domain(sequence, domain_format=spec.get("domain_format"), database="ALL")
    if domain["domain_count"] != 1:
        raise ValueError("Design requires one numbered domain per input sequence; split multidomain constructs explicitly")
    model = spec.get("profile")
    validate_profile(model)
    if model["species"] != "human":
        raise ValueError("Humanization/human-VH retention requires an explicitly human reference")
    parent = score_numbered(model, domain["numbering"], chain=domain["chain"], domain_format=domain["domain_format"])
    if parent["coverage"] != 1:
        raise ValueError("Incomplete fixed-model coordinate coverage")
    parent_provenance = spec.get("parent_provenance", {})
    if not isinstance(parent_provenance, dict):
        raise ValueError("parent_provenance must be an explicit JSON object")
    json.dumps(parent_provenance, allow_nan=False)
    if parent_provenance and parent_provenance.get("sequence_sha256") != domain["sequence_sha256"]:
        raise ValueError("Parent provenance must bind the exact input sequence SHA256")
    threshold = _threshold(spec["threshold"], model)
    protections = _protected(domain["numbering"], ())
    for pos in _vernier_positions(domain):
        protections.setdefault(pos, []).append("kabat_mapped_vernier")
    for field in ("protected_positions", "core_positions", "antigen_positions", "interface_positions"):
        for pos in _positions(spec.get(field, []), domain, field):
            if autonomy and field == "interface_positions":
                continue
            protections.setdefault(pos, []).append("caller_" + field)
    explicit = "allowed_substitutions" in spec
    if explicit == ("template" in spec):
        raise ValueError("Choose exactly one explicit allowed_substitutions map or one human V/J template")
    if explicit:
        alphabet, template_info = spec["allowed_substitutions"], None
    else:
        alphabet, template_info = template_alphabet(domain, spec["template"])
    if not isinstance(alphabet, dict):
        raise ValueError("allowed_substitutions must be a position-to-amino-acid-list map")
    mutable = None if "mutable_positions" not in spec else _positions(spec["mutable_positions"], domain, "mutable_positions")
    overrides = spec.get("protection_overrides", {})
    if overrides and not autonomy:
        raise ValueError("Humanization does not silently override conserved/Vernier protections")
    if not isinstance(overrides, dict):
        raise ValueError("protection_overrides must map explicit positions to review reasons")
    for pos, reason in overrides.items():
        if pos not in alphabet or not isinstance(reason, str) or not reason.strip():
            raise ValueError("Each autonomy override needs a declared allele position and a review reason")
        forbidden = [rule for rule in protections.get(pos, []) if rule not in {"kabat_mapped_vernier", "conservative_backbone_sensitive_residue"}]
        if forbidden:
            raise ValueError("Cannot override CDR, insertion, cysteine or caller immutable protection at " + pos)
    effective, excluded, relaxed = {}, {}, {}
    for pos, values in alphabet.items():
        _position(pos)
        if pos not in domain["numbering"] or not isinstance(values, (list, tuple)) or not values:
            raise ValueError("Every allele map position must exist and have a nonempty residue list")
        if any(not isinstance(aa, str) or aa not in AA for aa in values):
            raise ValueError("Replacement residues must be canonical uppercase amino acids")
        replacements = sorted(set(values) - {domain["numbering"][pos]})
        if not replacements:
            continue
        locks = protections.get(pos, [])
        if pos in overrides:
            relaxed[pos] = {"rules": list(locks), "reason": overrides[pos]}
            locks = []
        reasons = list(locks)
        if mutable is not None and pos not in mutable:
            reasons.append("outside_explicit_mutable_set")
        if reasons:
            if explicit:
                raise ValueError("Declared substitution conflicts with protection at %s: %s" % (pos, reasons))
            excluded[pos] = reasons
            continue
        if any(aa in "CGP" for aa in replacements) and (pos not in overrides or "C" in replacements):
            if not explicit:
                excluded[pos] = ["conservative_no_introduced_Cys_Gly_Pro"]
                continue
            raise ValueError("Introducing Cys is forbidden; Gly/Pro require an explicit autonomy override")
        effective[pos] = replacements
    guard_specs = []
    guards = spec.get("guards", [])
    if not isinstance(guards, (list, tuple)) or len(guards) > 8:
        raise ValueError("At most eight explicit guard profiles per domain")
    for guard in guards:
        if not isinstance(guard, dict) or set(guard) != {"profile", "threshold"}:
            raise ValueError("Each guard requires exactly profile and threshold")
        guard_model = guard["profile"]
        validate_profile(guard_model)
        score = score_numbered(guard_model, domain["numbering"], chain=domain["chain"], domain_format=domain["domain_format"])
        if score["coverage"] != 1:
            raise ValueError("Guard needs full matching chain/format coverage")
        guard_specs.append((guard_model, _threshold(guard["threshold"], guard_model)))
    return {"id": identifier, "domain": domain, "profile": model, "threshold": threshold, "guards": guard_specs,
            "parent_score": parent, "alphabet": effective, "protections": protections, "overrides": relaxed,
            "template": template_info, "excluded_template_positions": excluded,
            "baseline_liabilities": _liabilities(sequence), "parent_provenance": parent_provenance}


def _materialize(compiled, candidate):
    domain, identifier = compiled["domain"], compiled["id"]
    numbering = {pos: candidate[identifier + ":" + pos] for pos in domain["numbering"]}
    residues = list(domain["sequence"])
    for pos, index in domain["sequence_indices"].items():
        residues[index] = numbering[pos]
    return numbering, "".join(residues)


def _native_validator(compiled):
    def check(candidate):
        numbering, sequence = _materialize(compiled, candidate)
        result = _native_map(sequence)
        if result is None:
            return {"accepted": False, "evidence": {"reason": "successful_native_no_hit"}}
        original = compiled["domain"]
        preserved = all(result.get(key) == original[key] for key in ("sequence_indices", "chain", "domain_count", "start", "end")) and result.get("numbering") == numbering
        return {"accepted": preserved, "evidence": {"domain_id": compiled["id"], "numbering_preserved": preserved,
                                                     "native_result": result}}
    return check


def _scope(compiled, task, extra):
    from . import _native, __version__
    return {"task": task, "domains": [{"id": item["id"], "query": {key: value for key, value in item["domain"].items() if key != "sequence"},
                                         "parent_provenance": item["parent_provenance"],
                                         "profile_sha256": item["profile"]["model_sha256"], "threshold": item["threshold"],
                                         "guards": [{"profile_sha256": model["model_sha256"], "threshold": threshold} for model, threshold in item["guards"]],
                                         "protections": item["protections"], "overrides": item["overrides"],
                                         "template": item["template"], "excluded_template_positions": item["excluded_template_positions"],
                                         "baseline_motif_locations": [list(location) for location in sorted(item["baseline_liabilities"])]} for item in compiled],
            "numeric_policy": "binary64_fsum_mean_log_and_direct_midrank_inclusive_comparison_v1_no_epsilon",
            "motif_policy": "new_N_nonP_ST_N_GS_D_GS_occurrences_at_raw_indices_v1",
            "vernier_reference": "Foote-Winter1992;AbNumber-common-blob-d494cb3593745ee824727537176da5eec67e1aaf",
            "native_binary_sha256": hashlib.sha256(Path(_native.__file__).read_bytes()).hexdigest(),
            "compiler_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "acceptance_source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                         for name in ("antibody_design.py", "species_adaptation.py", "species_reference.py", "minimum_edit.py", "_api.py", "__init__.py")},
            "python_version": sys.version,
            "package_version": __version__, "numbering": {"database": "ALL", "scheme": "imgt", "ncpu": 1},
            "extra": extra}


def _run(compiled, *, task, max_mutations, extra_screen=None, extra_scope=None, **budgets):
    parent, allowed = {}, {}
    for item in compiled:
        parent.update({item["id"] + ":" + pos: aa for pos, aa in item["domain"]["numbering"].items()})
        allowed.update({item["id"] + ":" + pos: alleles for pos, alleles in item["alphabet"].items()})
    metadata = _scope(compiled, task, extra_scope or {})
    def screen(candidate):
        evidence = {"objectives": [], "new_motifs": {}}
        accepted = True
        for item in compiled:
            numbering, sequence = _materialize(item, candidate)
            ok, value = _objective(item["profile"], numbering, item["threshold"])
            evidence["objectives"].append({"domain_id": item["id"], "role": "target", **value})
            accepted = accepted and ok
            for model, threshold in item["guards"]:
                ok, value = _objective(model, numbering, threshold)
                evidence["objectives"].append({"domain_id": item["id"], "role": "guard", "profile_sha256": model["model_sha256"], **value})
                accepted = accepted and ok
            new = sorted(_liabilities(sequence) - item["baseline_liabilities"])
            evidence["new_motifs"][item["id"]] = [list(location) for location in new]
            accepted = accepted and not new
        if extra_screen is not None:
            ok, details = extra_screen(candidate)
            accepted = accepted and ok
            evidence["task_specific"] = details
        return {"accepted": bool(accepted), "evidence": evidence}
    validators = [({"kind": "native_exact_single_domain_mapping", "domain_id": item["id"],
                    "native_binary_sha256": metadata["native_binary_sha256"]}, _native_validator(item)) for item in compiled]
    result = solve_minimum(parent, allowed, max_edits=max_mutations, screen=screen, validators=validators,
                           scope_metadata=metadata, **budgets)
    result.update(task=task, solver_status=result["status"],
                  computational_screen_status="all_declared_checks_passed" if result["candidate"] is not None else "no_accepted_witness",
                  experimental_validation=False, evidence_quality="computational_hypothesis_only",
                  preflight_numbering_calls=2 * len(compiled), preflight_budget="bounded inputs; outside cooperative search deadline",
                  structural_core_assessment="not_inferred; caller protections only",
                  likelihood_family_conditioned=False, paired_statistical_model=False,
                  limitations=["Minimum is only for the declared edit alphabet, constraints and binary64 evaluator, not biological humanization/autonomy.",
                               "A cohort percentile is not species probability, immune tolerance or a validated humanization threshold.",
                               "CDR sequence preservation does not prove binding preservation; framework and partner effects remain unvalidated."])
    if result["candidate"] is not None:
        result["designed_domains"] = []
        for item in compiled:
            numbering, sequence = _materialize(item, result["candidate"])
            mutations = [{"position": pos, "scheme": "imgt", "sequence_index": item["domain"]["sequence_indices"][pos],
                          "from": aa, "to": numbering[pos]} for pos, aa in item["domain"]["numbering"].items() if numbering[pos] != aa]
            result["designed_domains"].append({"id": item["id"], "parent_sequence": item["domain"]["sequence"], "sequence": sequence,
                                               "mutations": mutations, "cdr_sequence_preserved": True,
                                               "parent_score": item["parent_score"],
                                               "score": score_numbered(item["profile"], numbering, chain=item["domain"]["chain"], domain_format=item["domain"]["domain_format"])})
    return result


def humanize(domains, *, max_mutations=5, **budgets):
    """Minimum-in-scope human reference target; one VH, one VL, or an explicit pair.

    Domain dicts contain id, sequence, domain_format, profile, threshold and exactly
    one allowed_substitutions map or explicit human V/J template. No model defaults.
    """
    domains = copy.deepcopy(domains)
    if not isinstance(domains, list) or not 1 <= len(domains) <= 2:
        raise ValueError("Supply one domain or an explicit VH/VL pair")
    if any(not isinstance(spec, dict) or not isinstance(spec.get("id"), str) for spec in domains):
        raise ValueError("Each domain must be an object with an explicit string id")
    if len({spec["id"] for spec in domains}) != len(domains):
        raise ValueError("Domain IDs must be unique")
    domains.sort(key=lambda spec: spec["id"])
    if any(spec.get("domain_format") not in {"VH", "VL"} for spec in domains):
        raise ValueError("Humanization mode supports conventional VH/VL, not an implicit VHH format bridge")
    if len(domains) == 2 and {spec["domain_format"] for spec in domains} != {"VH", "VL"}:
        raise ValueError("Paired mode requires one VH and one VL with separate compatible human profiles")
    compiled = [_compile_domain(spec) for spec in domains]
    return _run(compiled, task="minimum_edit_human_reference_target", max_mutations=max_mutations, **budgets)
