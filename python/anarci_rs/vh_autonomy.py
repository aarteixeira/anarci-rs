"""Explicit VH-interface engineering hypotheses; never a VH-to-VHH success label."""
import copy
import hashlib
import re
from pathlib import Path

from .antibody_design import HYDROPHOBIC, _compile_domain, _run, _materialize, _positions, _raw_score
from .species_reference import validate_profile
from .autonomy_structure import validate_interface_evidence


def project_vh_to_vhh(numbering, profile):
    """Descriptive cross-format projection with a fixed denominator, never a percentile."""
    validate_profile(profile)
    if profile["chain"] != "H" or profile["domain_format"] != "VHH":
        raise ValueError("A projection reference must explicitly be a VHH profile")
    missing = [pos for pos in profile["positions"] if pos not in numbering]
    invalid = [pos for pos in profile["positions"] if pos in numbering and numbering[pos] not in profile["probabilities"][pos]]
    if invalid:
        raise ValueError("Invalid residue in VHH projection")
    return {"query_format": "VH", "reference_format": "VHH", "kind": "cross_format_descriptive_proxy",
            "profile_id": profile["profile_id"], "model_sha256": profile["model_sha256"],
            "coverage": (len(profile["positions"]) - len(missing)) / len(profile["positions"]),
            "status": "incomplete_fixed_panel" if missing else "projected",
            "raw_score": None if missing else _raw_score(profile, numbering), "percentile": None,
            "missing_positions": missing, "reference_evidence": profile["evidence"],
            "meaning": "Not VHH membership, autonomy, solubility or biological fitness. Never an acceptance objective."}


def engineer_vh_autonomy(domain, *, hypothesis, hydrophobic_maximum,
                         interface_evidence=None, sequence_only=False,
                         interface_positions=None, interface_review=None,
                         projection_profiles=(), max_mutations=3, **budgets):
    """Minimum edits for a caller-declared, frozen interface hypothesis.

    The defining proxy is a ceiling on the count of AVILMFWY residues in the
    explicit interface set. This is NOT exposed surface area or an autonomy model.
    Optional VHH references are descriptive only, including if their score falls.
    """
    domain = copy.deepcopy(domain)
    if not isinstance(domain, dict) or domain.get("domain_format") != "VH":
        raise ValueError("Autonomy mode accepts conventional VH input, not a relabeled VHH")
    if "allowed_substitutions" not in domain or "template" in domain:
        raise ValueError("Autonomy hypotheses require an explicit reviewed substitution alphabet")
    if not isinstance(hypothesis, str) or not hypothesis.strip():
        raise ValueError("Describe the interface-engineering hypothesis")
    if type(sequence_only) is not bool:
        raise ValueError("sequence_only must be boolean")
    if interface_evidence is None and not sequence_only:
        return {"task": "vh_autonomy_hypothesis", "status": "needs_structure_evidence", "solver_status": "not_run",
                "computational_screen_status": "not_run", "evidence_quality": "insufficient_interface_evidence",
                "experimental_validation": False, "autonomous_function_established": False,
                "required": "Supply reviewed parent VH/VL evidence or explicitly select sequence_only with a reviewed interface annotation"}
    if not isinstance(interface_review, dict) or set(interface_review) != {"source_id", "review_note"} or any(not isinstance(v, str) or not v.strip() for v in interface_review.values()):
        raise ValueError("Explicit interface review requires source_id and review_note")
    if interface_evidence is not None and sequence_only:
        raise ValueError("Choose structural evidence or explicit sequence-only mode, not both")
    evidence = copy.deepcopy(interface_evidence)
    if evidence is not None:
        if (not isinstance(evidence, dict) or evidence.get("schema_version") != 1
                or evidence.get("query_format") != "VH"
                or evidence.get("structure_assessment") != "parent_CA_proximity_only"
                or evidence.get("interface_scope") != "frozen_parent"
                or evidence.get("provenance", {}).get("reviewed") is not True
                or not re.fullmatch(r"[0-9a-f]{64}", str(evidence.get("structure_sha256", "")))):
            raise ValueError("Unsupported or incomplete parent interface evidence")
        sequence = domain.get("sequence")
        if not isinstance(sequence, str) or evidence.get("query_sequence_sha256") != hashlib.sha256(sequence.encode()).hexdigest():
            raise ValueError("Interface evidence belongs to a different parent sequence")
        if interface_positions is not None:
            raise ValueError("Structural mode uses its frozen measured interface set; do not replace it")
        interface_positions = evidence["interface_positions"]
        domain["antigen_positions"] = sorted(set(domain.get("antigen_positions", [])) | set(evidence["antigen_protected_positions"]))
    if not isinstance(interface_positions, (list, tuple)) or not interface_positions or len(set(interface_positions)) != len(interface_positions):
        raise ValueError("Specify a nonempty unique frozen interface-position list")
    if type(hydrophobic_maximum) is not int or not 0 <= hydrophobic_maximum <= len(interface_positions):
        raise ValueError("hydrophobic_maximum is an explicit integer ceiling within the interface set size")
    domain["interface_positions"] = list(interface_positions)
    compiled = _compile_domain(domain, autonomy=True)
    query = compiled["domain"]
    if query["start"] != 0 or query["end"] != len(query["sequence"]) - 1:
        raise ValueError("Autonomy mode needs the isolated numbered VH sequence, not flanks or a full heavy-chain construct")
    interface = _positions(list(interface_positions), query, "interface")
    if not set(compiled["alphabet"]) <= interface:
        raise ValueError("Autonomy edit alphabet must be restricted to the reviewed interface set")
    if evidence is not None:
        expected = {key: query[key] for key in ("start", "end", "numbering", "sequence_indices")}
        if evidence.get("domain") != expected:
            raise ValueError("Interface evidence uses a different residue map")
        validate_interface_evidence(evidence, query["sequence"], expected)
    profiles = copy.deepcopy(list(projection_profiles))
    if len(profiles) > 8:
        raise ValueError("At most eight descriptive VHH projections")
    projections = [project_vh_to_vhh(query["numbering"], profile) for profile in profiles]
    if len({value["model_sha256"] for value in projections}) != len(projections):
        raise ValueError("Duplicate projection model")
    initial_count = sum(query["numbering"][pos] in HYDROPHOBIC for pos in interface)
    def screen(candidate):
        numbering, _ = _materialize(compiled, candidate)
        count = sum(numbering[pos] in HYDROPHOBIC for pos in interface)
        return count <= hydrophobic_maximum, {"metric": "frozen_interface_hydrophobic_residue_count",
                                              "alphabet": "".join(sorted(HYDROPHOBIC)), "parent": initial_count,
                                              "candidate": count, "maximum": hydrophobic_maximum,
                                              "is_surface_area": False, "is_solubility_prediction": False}
    scope = {"hypothesis": hypothesis, "interface_positions": sorted(interface), "interface_review": interface_review,
             "parent_interface_evidence": evidence, "sequence_only": sequence_only,
             "hydrophobic_alphabet": "".join(sorted(HYDROPHOBIC)), "hydrophobic_maximum": hydrophobic_maximum,
             "projection_profiles": [profile["model_sha256"] for profile in profiles],
             "autonomy_module_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "structure_evidence_module_sha256": hashlib.sha256(Path(__file__).with_name("autonomy_structure.py").read_bytes()).hexdigest()}
    result = _run([compiled], task="vh_autonomy_hypothesis", max_mutations=max_mutations,
                  extra_screen=screen, extra_scope=scope, **budgets)
    result["solver_status"] = result["status"]
    result.update(query_format="VH", output_format="engineered_VH_candidate", hypothesis=hypothesis,
                  autonomous_function_established=False, structural_validation=False,
                  evidence_mode="sequence_only_reviewed_annotation" if sequence_only else "parent_CA_proximity_only",
                  parent_vhh_projections=projections,
                  missing_evidence=["Candidate folding, solubility, stability, aggregation and binding are untested.",
                                    "Parent VL contributions to antigen binding may be lost."] +
                                   (["No parent structure or observed VH/VL interface was supplied."] if sequence_only else list(evidence.get("unknowns", []))),
                  experimental_next_steps=["Compare original paired Fv, isolated parent VH and variants in matched constructs/conditions.",
                                           "Measure expression/recovery, concentration-dependent monomer/aggregation, folding and stability.",
                                           "Measure affinity, specificity and the application-relevant functional assay; predefine acceptance criteria."])
    if result["candidate"] is not None:
        numbering, _ = _materialize(compiled, result["candidate"])
        candidate_projections = [project_vh_to_vhh(numbering, profile) for profile in profiles]
        for original, variant in zip(projections, candidate_projections):
            variant["raw_score_delta"] = None if original["raw_score"] is None or variant["raw_score"] is None else variant["raw_score"] - original["raw_score"]
        result["designed_domains"][0]["vhh_projections"] = candidate_projections
    return result
