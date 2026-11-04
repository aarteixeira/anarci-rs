"""Experimental, stdlib-only repertoire reference likelihoods (Python 3.9+).

These independent-position scores are not probabilities of species, biological
validation, or mutation recommendations. Provenance supplied by callers is not
independently verified. No reference data or trained defaults ship here.
"""

import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
FRAMEWORK_POSITIONS = tuple(str(i) for a, b in ((1, 26), (39, 55), (66, 104), (118, 128)) for i in range(a, b + 1))
_SPLITS = ("train", "calibration", "test")


def _json(value):
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise ValueError("metadata must be finite JSON data") from exc


def _hash(value):
    return hashlib.sha256(_json(value).encode("utf-8")).hexdigest()


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(name + " must be a nonempty string")
    return value


def _identity(species, chain, domain_format):
    _text(species, "species")
    if not re.fullmatch(r"[a-z][a-z0-9_]*", species):
        raise ValueError("species must be a canonical lowercase identifier")
    if chain not in ("H", "K", "L"):
        raise ValueError("chain must be H, K or L")
    if domain_format not in ("VH", "VL", "VHH"):
        raise ValueError("unknown domain_format")
    if (chain == "H") != (domain_format in ("VH", "VHH")):
        raise ValueError("incompatible chain and domain_format")


def _imgt_position_key(position):
    """Sequence order of ANARCI IMGT codes, independent of mapping order.

    IMGT inserts toward the centre from both anchors: right anchors 33, 61,
    and 112 run in reverse. Suffixes sort by length then alphabetically,
    which also preserves the native repeated-letter A..Z, AA..ZZ alphabet.
    """
    match = re.fullmatch(r"([1-9][0-9]{0,2})([A-Z]*)", position)
    if not match:
        raise ValueError("invalid IMGT position")
    base, suffix = int(match.group(1)), match.group(2)
    if not 1 <= base <= 128:
        raise ValueError("IMGT position outside 1..128")
    rank = 0
    for letter in suffix:
        rank = 26 * rank + ord(letter) - ord("A") + 1
    return base, -rank if base in (33, 61, 112) else rank


def _numbered_sequence(numbering):
    return "".join(numbering[p] for p in sorted(numbering, key=_imgt_position_key))


def _numbering(numbering):
    if not isinstance(numbering, dict):
        raise ValueError("numbering must be a position-to-residue dict")
    out = {}
    for position, residue in numbering.items():
        if not isinstance(position, str) or not re.fullmatch(r"([1-9][0-9]{0,2})([A-Z]*)", position):
            raise ValueError("invalid IMGT position")
        _imgt_position_key(position)
        if residue == "-":
            continue
        if not isinstance(residue, str) or len(residue) != 1 or residue not in AMINO_ACIDS:
            raise ValueError("numbering residues must be canonical amino acids or '-' gaps")
        out[position] = residue
    return out


def _positions(positions):
    if not isinstance(positions, (list, tuple)) or not positions:
        raise ValueError("positions must be a nonempty list of canonical IMGT positions")
    if any(not isinstance(p, str) or p not in {str(i) for i in range(1, 129)} for p in positions):
        raise ValueError("modeled positions must be canonical 1..128 without insertions")
    if len(set(positions)) != len(positions):
        raise ValueError("duplicate modeled positions")
    return sorted(positions, key=int)


def _groups(rows):
    return {
        "records": len(rows),
        "studies": len({r["study"] for r in rows}),
        "donors": len({(r["study"], r["donor"]) for r in rows}),
        "lineages": len({r["lineage"] for r in rows}),
    }


def _mean_score(probabilities, positions, numbering):
    terms = [math.log(probabilities[p][numbering[p]]) for p in positions if p in numbering]
    return math.fsum(terms) / len(terms) if terms else None


def _is_unresolved(value):
    return value.strip().lower().split(":", 1)[0] in {"unknown", "none", "no", "na", "n/a", "nan", "undefined", "missing", "null"}


def _unresolved_metadata(rows):
    return {key: sum(_is_unresolved(r[key]) for r in rows)
            for key in ("study", "donor", "lineage")}


def _donor_support(rows):
    donors = {(r["study"], r["donor"]) for r in rows}
    proxy = sum(_is_unresolved(donor) for _, donor in donors)
    return {"known_donors": len(donors) - proxy, "proxy_donors": proxy}


def _evidence(model):
    groups = model["audit"]["effective_groups"]
    units = model["source_manifest"].get("units", [])
    origin_declared = (isinstance(units, list) and bool(units) and all(
        isinstance(unit, dict) and unit.get("biological_origin") == "natural_repertoire"
        and unit.get("library_selection") == "unselected"
        and isinstance(unit.get("origin_evidence"), dict)
        and bool(unit["origin_evidence"].get("url"))
        and bool(unit["origin_evidence"].get("review_note")) for unit in units))
    supported = (origin_declared and len(model["positions"]) >= 60 and groups["train"]["studies"] >= 2
                 and groups["calibration"]["studies"] >= 1
                 and len(model["calibration_scores"]) >= 20 and bool(model["source_manifest"])
                 and not any(count for split in ("train", "calibration")
                             for count in model["audit"]["unresolved_metadata"][split].values()))
    return "repertoire_baseline" if supported else "research_only"


def fit_profile(records, *, species, chain, domain_format, atlas_id, version,
                positions=None, pseudocount=0.5, source_manifest=None):
    """Fit balanced smoothed frequencies; fail rather than repair split leakage.

    Study weight is equal; donors within study and lineages within donor receive
    equal weight; unique sequences within each lineage share its weight. Weights
    are scaled to the number of training studies before adding pseudocount per AA.
    Exact duplicate sequences within a split retain the lexicographically first
    metadata record (conflicting numbering fails). Gaps mean missing, never a 21st
    amino acid. Default positions require complete training coverage and >=60
    retained framework positions; explicit positions require complete coverage.
    Input sequence must equal the ungapped numbered domain in canonical IMGT
    order, with no flanks. Mapping hashes also enforce cross-split non-overlap.
    Calibration uses unweighted unique sequences from study-disjoint sources;
    these observations are not assumed independent within donors or lineages.
    """
    _identity(species, chain, domain_format)
    _text(atlas_id, "atlas_id")
    _text(version, "version")
    if isinstance(pseudocount, bool) or not isinstance(pseudocount, (int, float)) or not math.isfinite(pseudocount) or pseudocount <= 0:
        raise ValueError("pseudocount must be positive and finite")
    if source_manifest is None:
        source_manifest = {}
    if not isinstance(source_manifest, dict):
        raise ValueError("source_manifest must be a JSON object")
    source_manifest = json.loads(_json(source_manifest))
    if not isinstance(records, (list, tuple)) or not records:
        raise ValueError("records must be a nonempty list")
    rows = []
    owners = {key: {} for key in ("study", "donor", "lineage", "sequence", "numbering")}
    record_ids = set()
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("records must contain dicts")
        _json(record)
        for key in ("record_id", "species", "chain", "domain_format", "study", "donor", "lineage", "split", "sequence"):
            _text(record.get(key), key)
        if record["record_id"] in record_ids:
            raise ValueError("duplicate record_id")
        record_ids.add(record["record_id"])
        if (record["species"], record["chain"], record["domain_format"]) != (species, chain, domain_format):
            raise ValueError("record species/chain/domain_format mismatch")
        if record["split"] not in _SPLITS:
            raise ValueError("invalid split")
        if any(aa not in AMINO_ACIDS for aa in record["sequence"]):
            raise ValueError("sequence must contain only canonical amino acids")
        row = {k: record[k] for k in ("record_id", "study", "donor", "lineage", "split")}
        row["sequence_sha256"] = hashlib.sha256(record["sequence"].encode("ascii")).hexdigest()
        row["numbering"] = _numbering(record.get("numbering"))
        if record["sequence"] != _numbered_sequence(row["numbering"]):
            raise ValueError("sequence must exactly equal the reconstructed numbered domain (no flanks)")
        keys = {"study": row["study"], "donor": (row["study"], row["donor"]),
                "lineage": row["lineage"], "sequence": row["sequence_sha256"],
                "numbering": _hash(row["numbering"])}
        for kind, key in keys.items():
            previous = owners[kind].setdefault(key, row["split"])
            if previous != row["split"]:
                raise ValueError("split leakage: overlapping " + kind)
        rows.append(row)
    input_counts = dict(Counter(r["split"] for r in rows))
    unique = {}
    for row in sorted(rows, key=_json):
        key = row["sequence_sha256"]
        if key in unique and unique[key]["numbering"] != row["numbering"]:
            raise ValueError("duplicate sequence has conflicting numbering")
        unique.setdefault(key, row)
    rows = sorted(unique.values(), key=_json)
    by_split = {s: [r for r in rows if r["split"] == s] for s in _SPLITS}
    train = by_split["train"]
    if not train:
        raise ValueError("at least one training record required")
    requested = _positions(positions if positions is not None else FRAMEWORK_POSITIONS)
    selected = [p for p in requested if all(p in r["numbering"] for r in train)]
    if positions is not None and selected != requested:
        raise ValueError("explicit positions require complete training coverage")
    if not selected or (positions is None and len(selected) < 60):
        raise ValueError("default framework requires at least 60 fully covered positions")
    hierarchy = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for row in train:
        hierarchy[row["study"]][row["donor"]][row["lineage"]].append(row)
    weighted = []
    for study in sorted(hierarchy):
        donors = hierarchy[study]
        for donor in sorted(donors):
            lineages = donors[donor]
            for lineage in sorted(lineages):
                members = lineages[lineage]
                weight = 1.0 / (len(donors) * len(lineages) * len(members))
                weighted.extend((r, weight) for r in members)
    total = math.fsum(w for _, w in weighted)
    probabilities = {}
    for p in selected:
        probabilities[p] = {aa: (math.fsum(w for r, w in weighted if r["numbering"][p] == aa) + pseudocount)
                            / (total + 20 * pseudocount) for aa in AMINO_ACIDS}
    calibration = [r for r in by_split["calibration"] if all(p in r["numbering"] for p in selected)]
    scores = sorted(_mean_score(probabilities, selected, r["numbering"]) for r in calibration)
    provenance = [{k: v for k, v in r.items() if k != "numbering"} | {"numbering_sha256": _hash(r["numbering"])} for r in rows]
    model = {
        "schema_version": 1, "profile_id": species + ":" + chain + ":" + domain_format,
        "species": species, "chain": chain, "domain_format": domain_format,
        "atlas_id": atlas_id, "version": version, "positions": selected,
        "probabilities": probabilities, "calibration_scores": scores,
        "source_manifest": source_manifest, "record_provenance": provenance,
        "parameters": {"pseudocount": float(pseudocount), "requested_positions": requested,
                       "position_policy": "explicit" if positions is not None else "framework_intersection_min60",
                       "weighting": "equal_study_donor_lineage_unique_sequence", "gap_policy": "missing",
                       "calibration_weighting": "unweighted_unique_sequence"},
        "audit": {"input_counts": {s: input_counts.get(s, 0) for s in _SPLITS},
                  "exclusion_counts": {"exact_duplicates": len(records) - len(rows),
                                       "incomplete_calibration": len(by_split["calibration"]) - len(calibration)},
                  "effective_groups": {s: _groups(by_split[s]) for s in _SPLITS},
                  "unresolved_metadata": {s: _unresolved_metadata(by_split[s]) for s in _SPLITS},
                  "donor_support": {s: _donor_support(by_split[s]) for s in _SPLITS},
                  "full_coverage_calibration_groups": _groups(calibration),
                  "dropped_positions": [p for p in requested if p not in selected],
                  "training_weight_sum": total, "split_overlap": False},
    }
    model["evidence"] = _evidence(model)
    model["model_sha256"] = _hash(model)
    validate_profile(model)
    return model


def validate_profile(model):
    """Reject malformed or fingerprint-mismatched profiles (not authentication)."""
    try:
        if not isinstance(model, dict):
            raise ValueError("profile must be an object")
        _json(model)
        digest = model["model_sha256"]
        if not isinstance(digest, str) or digest != _hash({k: v for k, v in model.items() if k != "model_sha256"}):
            raise ValueError("model fingerprint mismatch")
        if model["schema_version"] != 1:
            raise ValueError("unsupported schema_version")
        _identity(model["species"], model["chain"], model["domain_format"])
        if model["profile_id"] != ":".join((model["species"], model["chain"], model["domain_format"])):
            raise ValueError("invalid profile_id")
        _text(model["atlas_id"], "atlas_id")
        _text(model["version"], "version")
        positions = _positions(model["positions"])
        if positions != model["positions"] or set(model["probabilities"]) != set(positions):
            raise ValueError("invalid position layout")
        for distribution in model["probabilities"].values():
            if not isinstance(distribution, dict) or set(distribution) != set(AMINO_ACIDS):
                raise ValueError("probabilities require exactly 20 amino acids")
            if any(isinstance(p, bool) or not isinstance(p, (int, float)) or not math.isfinite(p) or not 0 < p < 1 for p in distribution.values()):
                raise ValueError("invalid amino acid probability")
            if not math.isclose(math.fsum(distribution.values()), 1.0, rel_tol=0, abs_tol=1e-12):
                raise ValueError("probabilities must sum to one")
        scores = model["calibration_scores"]
        if not isinstance(scores, list) or any(isinstance(s, bool) or not isinstance(s, (int, float)) or not math.isfinite(s) or s >= 0 for s in scores):
            raise ValueError("invalid calibration scores")
        if scores != sorted(scores):
            raise ValueError("calibration scores must be sorted")
        if not isinstance(model["source_manifest"], dict) or not isinstance(model["record_provenance"], list):
            raise ValueError("invalid provenance")
        params = model["parameters"]
        if isinstance(params["pseudocount"], bool) or not isinstance(params["pseudocount"], (int, float)) or params["pseudocount"] <= 0:
            raise ValueError("invalid pseudocount")
        requested = _positions(params["requested_positions"])
        if not set(positions) <= set(requested):
            raise ValueError("positions not requested")
        if (params["gap_policy"] != "missing" or params["weighting"] != "equal_study_donor_lineage_unique_sequence"
                or params["calibration_weighting"] != "unweighted_unique_sequence"):
            raise ValueError("unsupported fit parameters")
        if params["position_policy"] not in ("explicit", "framework_intersection_min60"):
            raise ValueError("unsupported position policy")
        if params["position_policy"] == "framework_intersection_min60" and (len(positions) < 60 or requested != list(FRAMEWORK_POSITIONS)):
            raise ValueError("invalid default framework positions")
        if model["audit"]["split_overlap"] is not False:
            raise ValueError("split overlap is forbidden")
        audit = model["audit"]
        if set(audit["effective_groups"]) != set(_SPLITS) or set(audit["input_counts"]) != set(_SPLITS):
            raise ValueError("invalid split counts")
        groups_to_check = list(audit["effective_groups"].values()) + [audit["full_coverage_calibration_groups"]]
        for groups in groups_to_check:
            if set(groups) != {"records", "studies", "donors", "lineages"}:
                raise ValueError("invalid group layout")
            for count in groups.values():
                if type(count) is not int or count < 0:
                    raise ValueError("invalid group count")
            if not groups["studies"] <= groups["donors"] <= groups["records"] or groups["lineages"] > groups["records"]:
                raise ValueError("inconsistent group counts")
        for count in list(audit["input_counts"].values()) + list(audit["exclusion_counts"].values()):
            if type(count) is not int or count < 0:
                raise ValueError("invalid audit count")
        if audit["full_coverage_calibration_groups"]["records"] != len(scores):
            raise ValueError("inconsistent calibration count")
        provenance = model["record_provenance"]
        split_rows = {s: [] for s in _SPLITS}
        owners = {key: {} for key in ("study", "donor", "lineage", "sequence", "numbering")}
        record_ids = set()
        sequence_hashes = set()
        for row in provenance:
            for key in ("record_id", "study", "donor", "lineage"):
                _text(row[key], key)
            if row["split"] not in _SPLITS:
                raise ValueError("invalid provenance split")
            for key in ("sequence_sha256", "numbering_sha256"):
                if not isinstance(row[key], str) or not re.fullmatch(r"[0-9a-f]{64}", row[key]):
                    raise ValueError("invalid provenance hash")
            if row["record_id"] in record_ids or row["sequence_sha256"] in sequence_hashes:
                raise ValueError("duplicate provenance record")
            record_ids.add(row["record_id"])
            sequence_hashes.add(row["sequence_sha256"])
            keys = {"study": row["study"], "donor": (row["study"], row["donor"]),
                    "lineage": row["lineage"], "sequence": row["sequence_sha256"],
                    "numbering": row["numbering_sha256"]}
            for kind, key in keys.items():
                if owners[kind].setdefault(key, row["split"]) != row["split"]:
                    raise ValueError("provenance split leakage")
            split_rows[row["split"]].append(row)
        for split in _SPLITS:
            if _groups(split_rows[split]) != audit["effective_groups"][split] or audit["input_counts"][split] < len(split_rows[split]):
                raise ValueError("inconsistent provenance counts")
            if _unresolved_metadata(split_rows[split]) != audit["unresolved_metadata"][split]:
                raise ValueError("inconsistent unresolved metadata counts")
            if _donor_support(split_rows[split]) != audit["donor_support"][split]:
                raise ValueError("inconsistent donor support counts")
        if not split_rows["train"]:
            raise ValueError("missing training provenance")
        if sum(audit["input_counts"].values()) - len(provenance) != audit["exclusion_counts"]["exact_duplicates"]:
            raise ValueError("inconsistent duplicate count")
        if len(split_rows["calibration"]) - len(scores) != audit["exclusion_counts"]["incomplete_calibration"]:
            raise ValueError("inconsistent excluded calibration count")
        if audit["dropped_positions"] != [p for p in requested if p not in positions]:
            raise ValueError("inconsistent dropped positions")
        if not math.isclose(audit["training_weight_sum"], audit["effective_groups"]["train"]["studies"], rel_tol=0, abs_tol=1e-12):
            raise ValueError("invalid training weight sum")
        if model["evidence"] != _evidence(model):
            raise ValueError("invalid evidence tier")
    except (KeyError, TypeError, AttributeError, OverflowError) as exc:
        raise ValueError("malformed profile") from exc


def audit_profile(model):
    """Return counts and limits, never a biological-validation certificate."""
    validate_profile(model)
    result = json.loads(_json(model["audit"]))
    result.update({"profile_id": model["profile_id"], "model_sha256": model["model_sha256"],
                   "evidence": model["evidence"], "modeled_position_count": len(model["positions"]),
                   "reference_score_count": len(model["calibration_scores"]),
                   "provenance_verified": False, "biologically_validated": False})
    return result


def score_numbered(model, numbering, *, chain, domain_format):
    """Mean log likelihood and empirical midrank percentile on full coverage.

    Percentile is 0..100 relative to study-disjoint calibration sequences, not
    a species posterior. Unique sequences are unweighted, not assumed IID.
    No softmax or cross-profile normalization is applied.
    """
    validate_profile(model)
    _identity(model["species"], chain, domain_format)
    if (chain, domain_format) != (model["chain"], model["domain_format"]):
        raise ValueError("query chain/domain_format does not match profile")
    numbering = _numbering(numbering)
    contributions = []
    for p in model["positions"]:
        if p in numbering:
            aa = numbering[p]
            probability = model["probabilities"][p][aa]
            contributions.append({"position": p, "residue": aa, "probability": probability,
                                  "log_probability": math.log(probability)})
    coverage = len(contributions) / len(model["positions"])
    raw = _mean_score(model["probabilities"], model["positions"], numbering)
    scores = model["calibration_scores"]
    percentile = None
    warnings = ["Independent-position repertoire likelihood, not a probability of species or a validated mutation score.",
                "Caller-supplied source provenance has not been verified by this fitter.",
                "Calibration percentiles weight unique sequences equally; within-donor/lineage observations are not independent."]
    if coverage < 1:
        warnings.append("Partial modeled-position coverage: percentile unavailable.")
    elif not scores:
        warnings.append("No full-coverage study-disjoint calibration records: percentile unavailable.")
    else:
        percentile = 100.0 * (sum(s < raw for s in scores) + 0.5 * sum(s == raw for s in scores)) / len(scores)
    return {"status": "no_coverage" if coverage == 0 else "partial" if coverage < 1 else "ok",
            "raw_score": raw, "percentile": percentile, "coverage": coverage,
            "contributions": contributions, "warnings": warnings,
            "reference_score_count": len(scores), "support_status": model["evidence"],
            "calibration_weighting": model["parameters"]["calibration_weighting"],
            "calibration_groups": dict(model["audit"]["full_coverage_calibration_groups"]),
            "profile_id": model["profile_id"], "model_sha256": model["model_sha256"],
            "evidence": model["evidence"], "score_type": "mean_log_likelihood"}


def load_profile(path):
    """Load and validate a JSON profile without importing native extensions."""
    try:
        with Path(path).open(encoding="utf-8") as handle:
            model = json.load(handle)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError("invalid profile JSON") from exc
    validate_profile(model)
    return model


def save_profile(model, path):
    """Write deterministic JSON using exclusive creation; never overwrite."""
    validate_profile(model)
    with Path(path).open("x", encoding="utf-8") as handle:
        handle.write(_json(model) + "\n")
