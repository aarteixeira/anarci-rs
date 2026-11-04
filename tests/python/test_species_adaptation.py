"""Software invariants only; constructed profiles are not biological reference data.

Real OAS acquisition, fitting and held-out integration are separate from these
small numerical fixtures. Neither establishes biological adaptation success.
Native numbering tests use the existing public helper.
"""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from anarci_rs import species_adaptation as adaptation
from anarci_rs import species_cli
from anarci_rs.species_reference import fit_profile, save_profile
from safe_reference import HEAVY


def numerical_profile(species="human", residues=None, chain="H", training_studies=2):
    residues = residues or {"1": "Q", "2": "E", "23": "C", "27": "A", "41": "W"}
    records = []
    for index in range(8):
        split = "train" if index < 4 else "calibration" if index < 6 else "test"
        full_numbering = dict(residues, **{"128": "ACDEFGHI"[index]})
        domain_sequence = "".join(aa for pos, aa in sorted(full_numbering.items(), key=lambda item: int(item[0])))
        records.append({"record_id": "numeric-%d" % index, "species": species, "chain": chain, "domain_format": "VH" if chain == "H" else "VL",
                        "study": split + "-study-" + str(index % (training_studies if split == "train" else 2)), "donor": "donor-%d" % index,
                        "lineage": "numeric-lineage-%d" % index, "split": split,
                        "sequence": domain_sequence, "numbering": full_numbering})
    return fit_profile(records, species=species, chain=chain, domain_format="VH" if chain == "H" else "VL", atlas_id="software-test-only",
                       version="1", positions=list(residues))


class SpeciesAdaptationTests(unittest.TestCase):
    def test_independent_profiles_not_species_probabilities(self):
        query = {"1": "A", "2": "A", "23": "C", "27": "Y", "41": "W"}
        report = adaptation.profile_numbered(query, [numerical_profile(), numerical_profile("mouse")],
                                             chain="H", domain_format="VH", species=["human", "mouse", "dog"])
        self.assertEqual(report["unsupported_species"], ["dog"])
        self.assertEqual(len(report["profiles"]), 2)
        self.assertEqual(report["relative_preference"]["kind"], "matched_profile_log_score_ranking_not_probability")
        self.assertNotIn("species_probability", report["profiles"][0])
        self.assertEqual(report["profiles"][0]["evidence"], "research_only")

    def test_different_effective_smoothing_suppresses_relative_ranking(self):
        query = {"1": "A", "2": "A", "23": "C", "27": "Y", "41": "W"}
        result = adaptation.profile_numbered(query, [numerical_profile(), numerical_profile("mouse", training_studies=4)],
                                             chain="H", domain_format="VH")
        self.assertIsNone(result["relative_preference"])
        self.assertEqual(len(result["profiles"]), 2)

    def test_partial_queries_do_not_receive_calibrated_ranks(self):
        result = adaptation.profile_numbered({"1": "Q"}, [numerical_profile()], chain="H", domain_format="VH")
        self.assertIsNone(result["profiles"][0]["percentile"])
        self.assertLess(result["profiles"][0]["coverage"], 1)
        with self.assertRaises(ValueError):
            adaptation.propose_numbered({"1": "Q"}, numerical_profile(), chain="H", domain_format="VH")

    def test_no_silent_format_bridge(self):
        query = {"1": "A", "2": "A", "23": "C", "27": "Y", "41": "W"}
        result = adaptation.profile_numbered(query, [numerical_profile()], chain="H", domain_format="VHH")
        self.assertEqual(result["status"], "unsupported")
        with self.assertRaises(ValueError):
            adaptation.propose_numbered(query, numerical_profile(), chain="H", domain_format="VHH")
        with self.assertRaises(ValueError):
            adaptation.profile_numbered(query, [numerical_profile()], chain="B", domain_format="VH")

    def test_planner_protects_cdr_disulfide_and_caller_positions(self):
        query = {"1": "A", "2": "A", "23": "C", "27": "Y", "41": "W"}
        result = adaptation.propose_numbered(query, numerical_profile(), chain="H", domain_format="VH",
                                             protected_positions=(p for p in ["2"]), max_mutations=3)
        self.assertTrue(result["candidates"])
        for candidate in result["candidates"]:
            self.assertEqual([m["position"] for m in candidate["mutations"]], ["1"])
            self.assertGreater(candidate["raw_score_gain"], 0)
            self.assertFalse(candidate["experimental_validation"])
            self.assertEqual(candidate["evidence_level"], "score_improvement_only")
        self.assertIn("IMGT_CDR", result["protected_positions"]["27"])
        self.assertIn("caller_protected", result["protected_positions"]["2"])

    def test_explicit_mutable_set_and_budget_are_enforced(self):
        query = {"1": "A", "2": "A", "23": "C", "27": "Y", "41": "W"}
        result = adaptation.propose_numbered(query, numerical_profile(), chain="H", domain_format="VH",
                                             mutable_positions=["2"], max_evaluations=1)
        self.assertLessEqual(result["evaluations"], 1)
        self.assertEqual(result["candidates"][0]["mutations"][0]["position"], "2")
        self.assertEqual(adaptation.propose_numbered(query, numerical_profile(), chain="H", domain_format="VH",
                                                   mutable_positions=[])["status"], "no_eligible_improvement")
        for kwargs in ({"max_mutations": 0}, {"beam_width": True}, {"max_evaluations": 50001}, {"protected_positions": ["128"]}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                adaptation.propose_numbered(query, numerical_profile(), chain="H", domain_format="VH", **kwargs)

    def test_guard_scores_never_decrease(self):
        query = {"1": "A", "2": "A", "23": "C", "27": "Y", "41": "W"}
        guard = numerical_profile("mouse", query)
        result = adaptation.propose_numbered(query, numerical_profile(), chain="H", domain_format="VH", guard_profiles=[guard])
        self.assertEqual(result["status"], "no_eligible_improvement")

    def test_minimal_mutation_budgets_and_determinism(self):
        query = {"1": "A", "2": "A", "23": "C", "27": "Y", "41": "W"}
        first = adaptation.propose_numbered(query, numerical_profile(), chain="H", domain_format="VH")
        second = adaptation.propose_numbered(query, numerical_profile(), chain="H", domain_format="VH")
        self.assertEqual(first, second)
        self.assertEqual({c["mutation_count"] for c in first["candidates"]}, {1, 2})
        for row in first["candidates"]:
            self.assertEqual(row["mutation_count"], len(row["mutations"]))

    def test_native_position_mapping_preserves_public_sequence(self):
        domain = adaptation.number_domain(HEAVY, domain_format="VH")
        for position, residue in domain["numbering"].items():
            self.assertEqual(HEAVY[domain["sequence_indices"][position]], residue)
        self.assertEqual(domain["chain"], "H")
        self.assertEqual(domain["numbering_scheme"], "imgt")
        with self.assertRaises(ValueError):
            adaptation.number_domain(HEAVY, domain_format="VL")
        with self.assertRaises(ValueError):
            adaptation.number_domain("A" * 80, domain_format="VH")

    def test_mapping_rejects_native_coordinate_mismatch(self):
        with patch("anarci_rs.anarci", return_value=([[[[((1, " "), "Q")], 0, 0]]], [[{"chain_type": "H"}]], [])):
            with self.assertRaises(ValueError):
                adaptation.number_domain("A", domain_format="VH")

    def test_sequence_proposals_verify_raw_mutations_and_numbering(self):
        domain = adaptation.number_domain(HEAVY, domain_format="VH")
        pos = next(p for p, aa in domain["numbering"].items() if p in {"1", "2", "3"} and aa not in "CGP")
        original = domain["numbering"][pos]
        replacement = next(aa for aa in "QEAV" if aa != original)
        model = numerical_profile(residues={pos: replacement})
        result = adaptation.propose_sequence(HEAVY, model, domain_format="VH", mutable_positions=[pos])
        self.assertTrue(result["candidates"])
        for candidate in result["candidates"]:
            mutation = candidate["mutations"][0]
            self.assertEqual(HEAVY[mutation["sequence_index"]], mutation["from"])
            self.assertEqual(candidate["sequence"][mutation["sequence_index"]], mutation["to"])
            self.assertEqual(len(candidate["sequence"]), len(HEAVY))
            self.assertIsInstance(candidate["numbering_preserved"], bool)

    def test_light_chain_class_flip_is_not_eligible(self):
        parent = {"sequence": "A", "numbering": {"1": "A"}, "sequence_indices": {"1": 0},
                  "chain": "K", "domain_format": "VL", "domain_index": 0, "domain_count": 1}
        changed = dict(parent, chain="L", numbering={"1": "Q"})
        with patch.object(adaptation, "number_domain", side_effect=[parent, changed]):
            result = adaptation.propose_sequence("A", numerical_profile(residues={"1": "Q"}, chain="K"),
                                                 domain_format="VL", max_mutations=1)
        self.assertFalse(result["candidates"][0]["numbering_preserved"])
        self.assertFalse(result["candidates"][0]["eligible_for_review"])

    def test_cli_audit_and_output_no_clobber(self):
        with tempfile.TemporaryDirectory() as work:
            path = Path(work) / "model.json"
            output = Path(work) / "audit.json"
            save_profile(numerical_profile(), path)
            self.assertEqual(species_cli.main(["audit", str(path), "--output", str(output)]), 0)
            self.assertTrue(json.loads(output.read_text()))
            before = output.read_bytes()
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(species_cli.main(["audit", str(path), "--output", str(output)]), 2)
            self.assertEqual(before, output.read_bytes())


if __name__ == "__main__":
    unittest.main()
