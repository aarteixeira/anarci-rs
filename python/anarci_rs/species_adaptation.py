"""Research repertoire compatibility and conservative, score-only adaptation.

No reference model is downloaded or silently selected. Numbering remains native;
this transparent independent-position baseline does not predict biological fitness.
"""
import hashlib
import json
import math
import re

from .species_reference import score_numbered, validate_profile

AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
POSITION = re.compile(r"([1-9][0-9]{0,2})([A-Z]{0,4})\Z")
CDRS = ((27, 38), (56, 65), (105, 117))
LIMITS = [
    "Repertoire compatibility is not species identity, affinity, immunogenicity or organismal fitness.",
    "Independent-position profiles do not model epistasis, VH/VL pairing or three-dimensional structure.",
    "Percentiles are reference-cohort ranks, not a percentage of human/canine/camelid identity.",
    "Mutation proposals are computational hypotheses requiring independent evaluation and experiments.",
]


def _position(value):
    if not isinstance(value, str):
        raise ValueError("Positions must be IMGT strings, e.g. '42' or '111A'")
    match = POSITION.fullmatch(value)
    if not match or not 1 <= int(match[1]) <= 128:
        raise ValueError("Position must be within IMGT 1..128 with an optional insertion code")
    return int(match[1]), match[2]


def _numbering(value):
    if not isinstance(value, dict) or not value:
        raise ValueError("numbering must be a nonempty IMGT position-to-residue dictionary")
    result = {}
    for position, residue in value.items():
        _position(position)
        if not isinstance(residue, str) or residue not in AA:
            raise ValueError("Numbered residues must be canonical uppercase amino acids, without gaps")
        result[position] = residue
    return result


def profile_numbered(numbering, profiles, *, chain, domain_format, species=None):
    """Score each applicable species independently; do not manufacture class probabilities."""
    numbering = _numbering(numbering)
    if chain not in {"H", "K", "L"} or domain_format not in {"VH", "VL", "VHH"} or ((chain == "H") != (domain_format in {"VH", "VHH"})):
        raise ValueError("Specify a compatible antibody chain and domain format")
    profiles = list(profiles)
    if not profiles:
        raise ValueError("Provide explicit validated reference profiles")
    requested = None if species is None else list(species)
    if requested is not None and (not requested or any(not isinstance(s, str) or not s for s in requested)):
        raise ValueError("species must be a nonempty list of species identifiers")
    results, seen = [], set()
    compatible = []
    for model in profiles:
        validate_profile(model)
        if requested is not None and model["species"] not in requested:
            continue
        key = (model["species"], model["chain"], model["domain_format"])
        if key in seen:
            raise ValueError("Supply only one model version per species/chain/format in a report")
        seen.add(key)
        if model["chain"] != chain or model["domain_format"] != domain_format:
            continue
        score = score_numbered(model, numbering, chain=chain, domain_format=domain_format)
        results.append({"species": model["species"], "chain": chain, "domain_format": domain_format, **score})
        compatible.append(model)
    supported = {row["species"] for row in results}
    unsupported = [] if requested is None else sorted(set(requested) - supported)
    relative = None
    # Raw mean log scores can only be compared when the modeled coordinates and
    # estimator settings match. Even then the ranking is not a species posterior.
    if len(results) > 1 and all(row.get("coverage") == 1 for row in results):
        signatures = {(tuple(model["positions"]), json.dumps(model["parameters"], sort_keys=True), model["audit"]["training_weight_sum"]) for model in compatible}
        if len(signatures) == 1:
            relative = {"kind": "matched_profile_log_score_ranking_not_probability",
                        "ranking": [{"species": row["species"], "raw_score": row["raw_score"]}
                                    for row in sorted(results, key=lambda r: (-r["raw_score"], r["species"]))]}
    return {"schema_version": 1, "status": "scored" if results else "unsupported", "chain": chain,
            "domain_format": domain_format, "profiles": sorted(results, key=lambda r: r["species"]),
            "unsupported_species": unsupported, "relative_preference": relative,
            "warnings": list(LIMITS) + ([] if relative is not None or len(results) < 2 else ["Models have unmatched coordinates/settings or incomplete query coverage; no relative ranking reported."])}


def number_domain(sequence, *, domain_format, domain_index=0, database="ALL"):
    """Use the native engine and require an exact raw-sequence correspondence.

    domain_format is supplied by the caller: a lone H hit never establishes VHH.
    Coordinates in the returned map are zero-based input sequence indices.
    """
    if not isinstance(sequence, str) or not sequence or any(aa not in AA for aa in sequence):
        raise ValueError("sequence must be a nonempty uppercase canonical amino-acid string")
    if len(sequence) > 100_000:
        raise ValueError("sequence exceeds the native scan length limit")
    if domain_format not in {"VH", "VL", "VHH"}:
        raise ValueError("domain_format must be explicitly VH, VL or VHH")
    if type(domain_index) is not int or domain_index < 0:
        raise ValueError("domain_index must be a nonnegative integer")
    if database not in {"ALL", "pan"}:
        raise ValueError("This interface supports the pinned bundled ALL/pan numbering profiles")
    from . import anarci
    numbered, details, _ = anarci([("query", sequence)], scheme="imgt", database=database,
                                  allowed_species=None, ncpu=1)
    domains = numbered[0] or []
    if domain_index >= len(domains):
        raise ValueError("Requested numbered domain was not found; no species score can be assigned")
    positions, start, end = domains[domain_index]
    chain = details[0][domain_index]["chain_type"]
    if (domain_format in {"VH", "VHH"} and chain != "H") or (domain_format == "VL" and chain not in {"K", "L"}):
        raise ValueError("Requested domain format is incompatible with the numbered chain class")
    mapping, indices = {}, {}
    index = start
    for (position, insertion), residue in positions:
        if residue == "-":
            continue
        label = str(position) + insertion.strip()
        _position(label)
        if label in mapping or index is None or index >= len(sequence) or sequence[index] != residue:
            raise ValueError("Numbering does not provide a unique exact input-residue map")
        mapping[label], indices[label] = residue, index
        index += 1
    if not mapping or index != end + 1:
        raise ValueError("Numbering endpoints do not match the numbered residue sequence")
    return {"sequence": sequence, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "numbering": mapping, "sequence_indices": indices, "chain": chain, "domain_format": domain_format,
            "domain_index": domain_index, "domain_count": len(domains), "start": start, "end": end,
            "coordinate_convention": "zero-based input sequence; inclusive domain endpoints",
            "numbering_database": database, "numbering_scheme": "imgt",
            "format_provenance": "caller-supplied, not inferred from an unpaired heavy chain"}


def profile_sequence(sequence, profiles, *, domain_format, domain_index=0, database="ALL", species=None):
    domain = number_domain(sequence, domain_format=domain_format, domain_index=domain_index, database=database)
    result = profile_numbered(domain["numbering"], profiles, chain=domain["chain"], domain_format=domain_format, species=species)
    result["query"] = {key: value for key, value in domain.items() if key != "sequence"}
    return result


def _raw_score(model, numbering):
    # Internal hot path: callers have already validated each model and full query
    # coverage. All proposals are substitutions at existing canonical positions.
    positions = model["positions"]
    return math.fsum(math.log(model["probabilities"][pos][numbering[pos]]) for pos in positions) / len(positions)


def _bounded_integer(value, name, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError("%s must be an integer in %d..%d" % (name, low, high))


def _protected(numbering, extra):
    protections = {}
    for position, aa in numbering.items():
        base, insertion = _position(position)
        reasons = []
        if any(low <= base <= high for low, high in CDRS):
            reasons.append("IMGT_CDR")
        if insertion:
            reasons.append("insertion")
        if aa == "C" or base in {23, 104}:
            reasons.append("cysteine_or_canonical_disulfide_position")
        if aa in {"G", "P"}:
            reasons.append("conservative_backbone_sensitive_residue")
        if position in extra:
            reasons.append("caller_protected")
        if reasons:
            protections[position] = reasons
    return protections


def propose_numbered(numbering, target_profile, *, chain, domain_format,
                     protected_positions=(), mutable_positions=None, guard_profiles=(),
                     max_mutations=3, max_candidates=8, beam_width=16, max_evaluations=5000):
    """Generate bounded framework-only score-improving proposals, not validated therapeutics.

    Guard profiles are explicit same-format models whose raw score may not decrease.
    No cross-format VH/VHH bridge, CDR optimization or structural inference is implicit.
    """
    numbering = _numbering(numbering)
    validate_profile(target_profile)
    _bounded_integer(max_mutations, "max_mutations", 1, 5)
    _bounded_integer(max_candidates, "max_candidates", 1, 50)
    _bounded_integer(beam_width, "beam_width", 1, 64)
    _bounded_integer(max_evaluations, "max_evaluations", 1, 50_000)
    protected_positions = tuple(protected_positions)
    for position in protected_positions:
        _position(position)
        if position not in numbering:
            raise ValueError("Protected position is absent from the query: " + position)
    requested = None if mutable_positions is None else set(mutable_positions)
    if requested is not None:
        for position in requested:
            _position(position)
            if position not in numbering:
                raise ValueError("Mutable position is absent from the query: " + position)
    protections = _protected(numbering, set(protected_positions))
    parent = score_numbered(target_profile, numbering, chain=chain, domain_format=domain_format)
    if parent.get("coverage") != 1 or parent.get("raw_score") is None:
        raise ValueError("Adaptation requires complete coverage of the target model's fixed positions")
    guards = list(guard_profiles)
    guard_parent = []
    for model in guards:
        validate_profile(model)
        value = score_numbered(model, numbering, chain=chain, domain_format=domain_format)
        if value.get("coverage") != 1 or value.get("raw_score") is None:
            raise ValueError("Guard profiles must match the query chain/format and have complete coverage")
        guard_parent.append(value)
    options = []
    for position in target_profile["positions"]:
        if position in protections or (requested is not None and position not in requested):
            continue
        old = numbering[position]
        probabilities = target_profile["probabilities"][position]
        # Do not introduce Cys/Pro/Gly by default. This conservative heuristic is
        # not a substitute for structure, interface or Vernier-position review.
        for new in sorted(AA - {old, "C", "P", "G"}):
            delta = math.log(probabilities[new]) - math.log(probabilities[old])
            if delta > 1e-12:
                options.append((position, old, new, delta))
    options.sort(key=lambda item: (-item[3], _position(item[0]), item[2]))
    beam = [({}, numbering, parent)]
    candidates, visited, evaluations = [], set(), 0
    budget_exhausted = False
    for depth in range(1, max_mutations + 1):
        expanded = []
        for mutations, _, _ in beam:
            for position, old, new, local_delta in options:
                if position in mutations:
                    continue
                edits = dict(mutations, **{position: new})
                key = tuple(sorted(edits.items(), key=lambda item: _position(item[0])))
                if key in visited:
                    continue
                if evaluations >= max_evaluations:
                    budget_exhausted = True
                    break
                visited.add(key)
                variant = dict(numbering, **edits)
                score = {"raw_score": _raw_score(target_profile, variant)}
                evaluations += 1
                guard_scores = [{"raw_score": _raw_score(model, variant)} for model in guards]
                if any(value["raw_score"] + 1e-12 < initial["raw_score"] for value, initial in zip(guard_scores, guard_parent)):
                    continue
                gain = score["raw_score"] - parent["raw_score"]
                if gain <= 1e-12:
                    continue
                row = {"numbering": variant, "mutation_count": depth,
                       "mutations": [{"position": pos, "scheme": "imgt", "from": numbering[pos], "to": aa,
                                      "log_probability_gain": math.log(target_profile["probabilities"][pos][aa]) - math.log(target_profile["probabilities"][pos][numbering[pos]])}
                                     for pos, aa in key],
                       "score": score, "raw_score_gain": gain, "guard_scores": guard_scores,
                       "evidence_level": "score_improvement_only", "experimental_validation": False}
                candidates.append(row)
                expanded.append((edits, variant, score))
            if budget_exhausted:
                break
        expanded.sort(key=lambda item: (-item[2]["raw_score"], tuple(sorted(item[0].items()))))
        beam = expanded[:beam_width]
        if not beam or budget_exhausted:
            break
    # Preserve a solution at each mutation budget before filling with alternatives.
    chosen = []
    candidates.sort(key=lambda row: (row["mutation_count"], -row["raw_score_gain"], str(row["mutations"])))
    for depth in range(1, max_mutations + 1):
        best = next((row for row in candidates if row["mutation_count"] == depth), None)
        if best is not None and len(chosen) < max_candidates:
            chosen.append(best)
    for row in sorted(candidates, key=lambda value: (-value["raw_score_gain"], value["mutation_count"], str(value["mutations"]))):
        if len(chosen) >= max_candidates:
            break
        if row not in chosen:
            chosen.append(row)
    for rank, row in enumerate(chosen, 1):
        row["candidate_id"] = "candidate-%03d" % rank
        row["score"] = score_numbered(target_profile, row["numbering"], chain=chain, domain_format=domain_format)
        row["guard_scores"] = [score_numbered(model, row["numbering"], chain=chain, domain_format=domain_format) for model in guards]
    return {"schema_version": 1, "status": "proposed" if chosen else "no_eligible_improvement",
            "target_species": target_profile["species"], "chain": chain, "domain_format": domain_format,
            "parent_score": parent, "candidates": chosen, "protected_positions": protections,
            "evaluations": evaluations, "budget_exhausted": budget_exhausted,
            "search": {"method": "bounded_beam_search_independent_position_baseline", "max_mutations": max_mutations,
                       "beam_width": beam_width, "max_candidates": max_candidates, "max_evaluations": max_evaluations,
                       "globally_optimal": False, "epistasis_modeled": False},
            "warnings": list(LIMITS) + ["Framework surface exposure, Vernier roles, antigen contacts and VH/VL interfaces are not inferred. Supply protected/mutable positions from independent evidence.",
                                        "Preserving CDR sequence does not establish binding retention. Guard scores are not independent functional validation."]}


def _liabilities(sequence):
    motifs = {"potential_N_glycosylation": r"N[^P][ST]", "deamidation_prone_motif": r"N[GS]", "isomerization_prone_motif": r"D[GS]"}
    return {(name, match.start()) for name, pattern in motifs.items()
            for match in re.finditer("(?=(" + pattern + "))", sequence)}


def propose_sequence(sequence, target_profile, *, domain_format, domain_index=0, database="ALL", **options):
    """Number a sequence, propose changes, and map each substitution to the raw input."""
    domain = number_domain(sequence, domain_format=domain_format, domain_index=domain_index, database=database)
    result = propose_numbered(domain["numbering"], target_profile, chain=domain["chain"], domain_format=domain_format, **options)
    baseline_liabilities = _liabilities(sequence)
    for candidate in result["candidates"]:
        residues = list(sequence)
        for mutation in candidate["mutations"]:
            index = domain["sequence_indices"][mutation["position"]]
            if residues[index] != mutation["from"]:
                raise ValueError("Mutation wild-type residue does not match the input")
            residues[index] = mutation["to"]
            mutation["sequence_index"] = index
        variant = "".join(residues)
        candidate["sequence"] = variant
        candidate["sequence_sha256"] = hashlib.sha256(variant.encode()).hexdigest()
        candidate["new_liability_flags"] = [{"kind": kind, "sequence_index": index} for kind, index in sorted(_liabilities(variant) - baseline_liabilities)]
        # Renumber to ensure substitutions did not silently change domain boundaries.
        try:
            mapped = number_domain(variant, domain_format=domain_format, domain_index=domain_index, database=database)
            candidate["numbering_preserved"] = (mapped["sequence_indices"] == domain["sequence_indices"]
                                                and mapped["domain_count"] == domain["domain_count"]
                                                and mapped["chain"] == domain["chain"]
                                                and mapped["numbering"] == candidate["numbering"])
        except (ValueError, AssertionError):
            candidate["numbering_preserved"] = False
        candidate["eligible_for_review"] = candidate["numbering_preserved"] and not candidate["new_liability_flags"]
    result["query"] = {key: value for key, value in domain.items() if key != "sequence"}
    return result
