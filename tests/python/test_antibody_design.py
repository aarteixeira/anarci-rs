"""Algorithm/interface fixtures only, not evidence of humanization or autonomy."""
import copy
import hashlib
import math
import unittest
from unittest.mock import patch

from anarci_rs import antibody_design as design
from anarci_rs import vh_autonomy as autonomy
from anarci_rs.minimum_edit import EvaluationError
from anarci_rs.species_reference import fit_profile


def profile(values, chain="H", domain_format="VH", species="human"):
    rows = []
    for index in range(8):
        split = "train" if index < 4 else "calibration" if index < 6 else "test"
        numbering = dict(values, **{"128": "ACDEFGHI"[index]})
        rows.append({"record_id": str(index), "species": species, "chain": chain, "domain_format": domain_format,
                     "study": split + str(index % 2), "donor": str(index), "lineage": str(index), "split": split,
                     "sequence": "".join(numbering[pos] for pos in sorted(numbering, key=int)), "numbering": numbering})
    return fit_profile(rows, species=species, chain=chain, domain_format=domain_format, atlas_id="numeric-test-only", version="1", positions=list(values))


def domain(sequence, *, domain_format, **kwargs):
    return {"sequence": sequence, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "numbering": {str(i + 1): aa for i, aa in enumerate(sequence)},
            "sequence_indices": {str(i + 1): i for i in range(len(sequence))},
            "chain": "H" if domain_format == "VH" else "K", "domain_format": domain_format,
            "domain_count": 1, "domain_index": 0, "start": 0, "end": len(sequence) - 1}


def native(sequence, scheme="imgt"):
    # Artificial index correspondence tests logic, not a biological Kabat alignment.
    mapped = domain(sequence, domain_format="VH")
    if scheme == "kabat":
        mapped["numbering"] = {str(i + 5): aa for i, aa in enumerate(sequence)}
        mapped["sequence_indices"] = {str(i + 5): i for i in range(len(sequence))}
    return mapped


def request(sequence="AA", values=None, allowed=None):
    model = profile(values or {"1": "Q", "2": "E"})
    return {"id": "H", "sequence": sequence, "domain_format": "VH", "profile": model,
            "threshold": {"kind": "raw_score", "minimum": design._raw_score(model, values or {"1": "Q", "2": "E"}), "rationale": "numerical exact-target regression, not a biological cutoff"},
            "allowed_substitutions": allowed if allowed is not None else {"1": ["Q"], "2": ["E"]}}


class DesignTests(unittest.TestCase):
    def setUp(self):
        self.number = patch.object(design, "number_domain", side_effect=domain).start()
        self.native = patch.object(design, "_native_map", side_effect=native).start()
        self.addCleanup(patch.stopall)

    def test_true_minimum_with_native_validation_and_zero_edit(self):
        spec = request()
        result = design.humanize([spec], max_mutations=2)
        self.assertEqual(result["status"], "proven_minimum")
        self.assertEqual(result["minimum_edits"], 2)
        self.assertEqual(result["completed_excluded_through"], 1)
        self.assertEqual(result["designed_domains"][0]["sequence"], "QE")
        self.assertFalse(result["experimental_validation"])
        spec["sequence"] = "QE"
        result = design.humanize([spec], max_mutations=0)
        self.assertEqual(result["status"], "parent_meets_target")
        self.assertEqual(result["minimum_edits"], 0)

    def test_native_no_hit_rejects_but_exception_is_error(self):
        def no_hit(sequence, scheme="imgt"):
            return native(sequence, scheme) if scheme == "kabat" else None
        with patch.object(design, "_native_map", side_effect=no_hit):
            result = design.humanize([request()], max_mutations=2)
        self.assertEqual(result["status"], "proven_infeasible_within_scope")
        def broken(sequence, scheme="imgt"):
            if scheme == "kabat": return native(sequence, scheme)
            raise AssertionError("backend operational fault")
        with patch.object(design, "_native_map", side_effect=broken), self.assertRaises(EvaluationError):
            design.humanize([request()], max_mutations=2)

    def test_final_renumbering_is_inside_acceptance_not_postfilter(self):
        def remap(sequence, scheme="imgt"):
            result = native(sequence, scheme)
            if scheme == "imgt": result["chain"] = "K"
            return result
        with patch.object(design, "_native_map", side_effect=remap):
            result = design.humanize([request()], max_mutations=2)
        self.assertEqual(result["status"], "proven_infeasible_within_scope")

    def test_percentile_midrank_and_adjacent_threshold(self):
        spec = request()
        spec["threshold"] = {"kind": "percentile", "minimum": 50.0, "rationale": "tie test"}
        self.assertEqual(design.humanize([spec], max_mutations=2)["minimum_edits"], 2)
        spec["threshold"]["minimum"] = math.nextafter(50.0, math.inf)
        self.assertEqual(design.humanize([spec], max_mutations=2)["status"], "proven_infeasible_within_scope")

    def test_budget_exhaustion_never_claims_infeasible(self):
        result = design.humanize([request()], max_mutations=2, max_evaluations=1)
        self.assertEqual(result["status"], "budget_exhausted")
        self.assertFalse(result["proof"])
        self.assertIsNone(result["candidate"])

    def test_protection_conflicts_unknown_options_and_invalid_thresholds(self):
        for extra in ({"core_positions": ["1"]}, {"interface_positions": ["1"]}, {"antigen_positions": ["1"]},
                      {"protected_positions": ["1"]}, {"protected_position": ["1"]},
                      {"protection_overrides": {"1": "not allowed in A"}}):
            with self.subTest(extra=extra), self.assertRaises(ValueError):
                design.humanize([{**request(), **extra}])
        for value in (True, math.nan, math.inf):
            spec = request(); spec["threshold"]["minimum"] = value
            with self.assertRaises(ValueError): design.humanize([spec])

    def test_pair_minimizes_total_changes_and_requires_light_reference(self):
        h = request("A", {"1": "Q"}, {"1": ["Q"]})
        l = request("A", {"1": "E"}, {"1": ["E"]})
        l.update(id="L", domain_format="VL", profile=profile({"1": "E"}, chain="K", domain_format="VL"))
        l["threshold"]["minimum"] = design._raw_score(l["profile"], {"1": "E"})
        with patch.object(design, "_vernier_positions", return_value=set()), patch.object(design, "_native_validator", return_value=lambda c: {"accepted": True, "evidence": {"synthetic_validator": True}}):
            result = design.humanize([h, l], max_mutations=2)
            self.assertEqual(result["minimum_edits"], 2)
            self.assertEqual([d["sequence"] for d in result["designed_domains"]], ["Q", "E"])
            l["profile"] = h["profile"]
            with self.assertRaises(ValueError): design.humanize([h, l])

    def test_vernier_is_mapped_via_raw_index_not_imgt_number(self):
        def kabat(sequence, scheme="imgt"):
            result = native(sequence, scheme)
            if scheme == "kabat": result.update(numbering={"2": sequence[0], "5": sequence[1]}, sequence_indices={"2": 0, "5": 1})
            return result
        with patch.object(design, "_native_map", side_effect=kabat), self.assertRaisesRegex(ValueError, "vernier"):
            design.humanize([request()])

    def test_actual_template_contract_no_implicit_nearest_or_mosaic(self):
        from anarci_rs.germlines import all_germlines
        query = domain("AA", domain_format="VH")
        template = {"v_gene": "IGHV3-23*01", "j_gene": "IGHJ4*02"}
        alphabet, evidence = design.template_alphabet(query, template)
        self.assertEqual(evidence["semantics"], "single_V_J_template_substitutions")
        self.assertFalse(evidence["likelihood_family_conditioned"])
        self.assertEqual(alphabet.get("1"), [all_germlines["V"]["H"]["human"][template["v_gene"]][0]])
        with self.assertRaises(ValueError): design.template_alphabet(query, {"v_gene": "nearest"})

    def test_scope_and_run_hashes_have_distinct_meanings(self):
        a = design.humanize([request()], max_mutations=2, max_nodes=100)
        b = design.humanize([request()], max_mutations=2, max_nodes=200)
        self.assertEqual(a["scope_sha256"], b["scope_sha256"])
        self.assertNotEqual(a["run_sha256"], b["run_sha256"])
        spec = request(); spec["threshold"]["minimum"] -= 0.01
        self.assertNotEqual(a["scope_sha256"], design.humanize([spec], max_mutations=2)["scope_sha256"])


class NativeResultTests(unittest.TestCase):
    def test_malformed_raw_native_output_is_never_a_scientific_rejection(self):
        valid_domain = ([((1, " "), "A")], 0, 0)
        for numbered, details in (([[]], [None]), (["xx"], [None]),
                                  ([None], [[{"chain_type": "H"}]]),
                                  ([[valid_domain]], [[{"chain_type": None}]]),
                                  ([[valid_domain, valid_domain]], [[{"chain_type": "H"}]]),
                                  ([[valid_domain]], ["not-details"])):
            with self.subTest(numbered=numbered, details=details), patch("anarci_rs.anarci", return_value=(numbered, details, [[]])):
                with self.assertRaises(ValueError): design._native_map("A")
        with patch("anarci_rs.anarci", return_value=([None], [None], [[]])):
            self.assertIsNone(design._native_map("A"))
        with patch("anarci_rs.anarci", return_value=([[valid_domain, ([((1, " "), "A")], 1, 1)]], [[{"chain_type": "H"}, {"chain_type": "H"}]], [[]])):
            self.assertEqual(design._native_map("AA"), {"domain_count": 2})


class AutonomyTests(unittest.TestCase):
    setUp = DesignTests.setUp

    def arguments(self):
        spec = request("VA", {"1": "Q", "2": "A"}, {"1": ["Q"]})
        spec["threshold"]["minimum"] = -100.0
        return dict(domain=spec, hypothesis="Numerical hypothesis: reduce hydrophobic count on a declared site, not a biological recipe",
                    hydrophobic_maximum=0, sequence_only=True, interface_positions=["1"],
                    interface_review={"source_id": "synthetic-test-only", "review_note": "algorithm contract fixture"})

    def test_autonomy_is_an_explicit_hypothesis_with_real_solver(self):
        result = autonomy.engineer_vh_autonomy(**self.arguments())
        self.assertEqual(result["status"], "proven_minimum")
        self.assertEqual(result["minimum_edits"], 1)
        self.assertEqual(result["designed_domains"][0]["sequence"], "QA")
        self.assertFalse(result["autonomous_function_established"])
        self.assertFalse(result["structural_validation"])
        self.assertEqual(result["query_format"], "VH")
        self.assertEqual(result["evidence_mode"], "sequence_only_reviewed_annotation")

    def test_vhh_projection_can_decrease_and_never_provides_percentile(self):
        args = self.arguments()
        args["projection_profiles"] = [profile({"1": "V", "2": "A"}, domain_format="VHH", species="alpaca")]
        result = autonomy.engineer_vh_autonomy(**args)
        projection = result["designed_domains"][0]["vhh_projections"][0]
        self.assertLess(projection["raw_score_delta"], 0)
        self.assertIsNone(projection["percentile"])
        self.assertEqual(projection["query_format"], "VH")
        partial = autonomy.project_vh_to_vhh({"1": "V"}, args["projection_profiles"][0])
        self.assertIsNone(partial["raw_score"])
        self.assertEqual(partial["status"], "incomplete_fixed_panel")

    def test_missing_structure_not_converted_to_measured_interface(self):
        args = self.arguments(); args["sequence_only"] = False
        result = autonomy.engineer_vh_autonomy(**args)
        self.assertEqual(result["status"], "needs_structure_evidence")
        self.assertEqual(result["solver_status"], "not_run")

    def test_evidence_sequence_and_coordinate_mismatch_rejected(self):
        args = self.arguments(); args["sequence_only"] = False; args.pop("interface_positions")
        args["interface_evidence"] = {"schema_version": 1, "query_format": "VH", "structure_assessment": "parent_CA_proximity_only",
                                      "interface_scope": "frozen_parent", "provenance": {"reviewed": True},
                                      "structure_sha256": "a" * 64, "query_sequence_sha256": "b" * 64}
        with self.assertRaisesRegex(ValueError, "different parent"):
            autonomy.engineer_vh_autonomy(**args)

    def test_mutation_outside_interface_and_protection_override_fail_closed(self):
        args = self.arguments(); args["interface_positions"] = ["2"]
        with self.assertRaisesRegex(ValueError, "interface set"):
            autonomy.engineer_vh_autonomy(**args)
        args = self.arguments(); args["domain"]["antigen_positions"] = ["1"]
        args["domain"]["protection_overrides"] = {"1": "cannot bypass antigen locks"}
        with self.assertRaises(ValueError): autonomy.engineer_vh_autonomy(**args)


if __name__ == "__main__":
    unittest.main()
