"""Synthetic numerical/software invariants ONLY; not species evidence."""
import copy
import importlib.util
import math
from pathlib import Path
import tempfile
import unittest

# Standalone import deliberately works without a compiled native extension.
_SPEC = importlib.util.spec_from_file_location(
    "species_reference_under_test", Path(__file__).resolve().parents[2] / "python/anarci_rs/species_reference.py")
sr = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(sr)


def record(identifier, residue="A", split="train", study=None, donor="d", lineage=None):
    # Distinct unmodeled residues prevent exact-domain leakage in tiny fixtures.
    tag = sr.AMINO_ACIDS[len(identifier)]
    return {"record_id": identifier, "species": "human", "chain": "H", "domain_format": "VH",
            "study": study or split, "donor": donor, "lineage": lineage or identifier,
            "split": split, "sequence": residue + "C" + tag,
            "numbering": {"1": residue, "2": "C", "128": tag}}


def fit(rows, **kwargs):
    options = dict(species="human", chain="H", domain_format="VH", atlas_id="synthetic-invariants",
                   version="1", positions=["1", "2"])
    options.update(kwargs)
    return sr.fit_profile(rows, **options)


class SpeciesReferenceTests(unittest.TestCase):
    def test_arithmetic_and_midrank(self):
        model = fit([record("a"), record("bb", "C"), record("ccc", "A", "calibration"),
                     record("dddd", "C", "calibration")])
        self.assertAlmostEqual(model["probabilities"]["1"]["A"], 1 / 11)
        self.assertAlmostEqual(model["probabilities"]["2"]["C"], 1.5 / 11)
        score = sr.score_numbered(model, {"1": "A", "2": "C"}, chain="H", domain_format="VH")
        self.assertAlmostEqual(score["raw_score"], (math.log(1 / 11) + math.log(1.5 / 11)) / 2)
        self.assertEqual(score["percentile"], 50)
        self.assertEqual(score["evidence"], "research_only")
        self.assertEqual(len(score["contributions"]), 2)

    def test_partial_and_empty_queries(self):
        model = fit([record("a"), record("bb", split="calibration")])
        for numbering, status, coverage in [({"1": "A", "2": "-"}, "partial", .5), ({}, "no_coverage", 0)]:
            result = sr.score_numbered(model, numbering, chain="H", domain_format="VH")
            self.assertEqual((result["status"], result["coverage"]), (status, coverage))
            self.assertIsNone(result["percentile"])
        self.assertIsNone(sr.score_numbered(model, {}, chain="H", domain_format="VH")["raw_score"])

    def test_duplicate_and_order_independence(self):
        a, b = record("a"), record("bb", "C")
        duplicate = dict(a, record_id="other", donor="other", lineage="other")
        model = fit([a, b, duplicate])
        self.assertEqual(model, fit([duplicate, b, a]))
        self.assertEqual(model["probabilities"], fit([a, b])["probabilities"])
        self.assertEqual(model["audit"]["exclusion_counts"]["exact_duplicates"], 1)
        self.assertNotIn('"sequence":', sr._json(model))

    def test_balance_donors_lineages_and_studies(self):
        rows = [record("a", "A", donor="one"), record("bb", "C", donor="two"),
                record("ccc", "C", donor="two", lineage="bb")]
        p = fit(rows)["probabilities"]["1"]
        self.assertEqual(p["A"], p["C"])
        rows.append(record("dddd", "D", study="second"))
        p = fit(rows)["probabilities"]["1"]
        self.assertAlmostEqual(p["D"], 1.5 / 12)
        self.assertAlmostEqual(p["A"], 1 / 12)

    def test_leakage_rejected_before_dedup(self):
        base = record("a")
        for field in ("study", "lineage", "sequence"):
            other = record("bb", "C", "calibration")
            other[field] = base[field]
            if field == "sequence":
                other["numbering"] = dict(base["numbering"])
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "split leakage"):
                fit([base, other])

    def test_qualified_donors_permitted(self):
        model = fit([record("a", study="first"), record("bb", "C", "calibration", study="second")])
        self.assertEqual(model["audit"]["effective_groups"]["calibration"]["donors"], 1)

    def test_invalid_metadata_and_residues(self):
        for field, value in [("species", "mouse"), ("chain", "K"), ("domain_format", "VHH"),
                             ("study", float("nan")), ("split", "valid"), ("sequence", "ABX"),
                             ("numbering", {"1": "AA"}), ("numbering", {"129": "A"}),
                             ("numbering", {"01": "A"}), ("numbering", {"1": "X"})]:
            row = record("a")
            row[field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                fit([row])
        for kwargs in ({"chain": "T"}, {"domain_format": "unknown"}, {"chain": "K", "domain_format": "VH"},
                       {"pseudocount": float("inf")}, {"pseudocount": 0}, {"source_manifest": []}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                fit([record("a")], **kwargs)

    def test_query_validation(self):
        model = fit([record("a")])
        for kwargs in ({"chain": "K", "domain_format": "VL"}, {"chain": "H", "domain_format": "VHH"},
                       {"chain": "X", "domain_format": "VH"}):
            with self.assertRaises(ValueError):
                sr.score_numbered(model, {"1": "A"}, **kwargs)
        with self.assertRaises(ValueError):
            sr.score_numbered(model, {"1": "AC"}, chain="H", domain_format="VH")

    def test_default_framework_and_gap_policy(self):
        row = record("a")
        row["numbering"] = {p: "A" for p in sr.FRAMEWORK_POSITIONS}
        row["numbering"].update({"1": "-", "111A": "C"})
        row["sequence"] = "A" * 81 + "C" + "A" * 11
        model = fit([row], positions=None)
        self.assertEqual(len(model["positions"]), 92)
        self.assertEqual(model["audit"]["dropped_positions"], ["1"])
        self.assertNotIn("111A", model["positions"])
        with self.assertRaises(ValueError):
            fit([record("a")], positions=None)
        with self.assertRaises(ValueError):
            fit([row], positions=["1"])
        with self.assertRaises(ValueError):
            fit([row], positions=["111A"])

    def test_incomplete_calibration_excluded(self):
        row = record("bb", split="calibration")
        row["numbering"] = {"1": "A"}
        row["sequence"] = "A"
        model = fit([record("a"), row])
        self.assertEqual(model["calibration_scores"], [])
        self.assertEqual(sr.audit_profile(model)["exclusion_counts"]["incomplete_calibration"], 1)

    def test_tampering_and_malformed_probabilities(self):
        model = fit([record("a")])
        modified = copy.deepcopy(model)
        modified["probabilities"]["1"]["A"] = .2
        with self.assertRaisesRegex(ValueError, "fingerprint"):
            sr.validate_profile(modified)
        for distribution in ({"A": 1.0}, {aa: .1 for aa in sr.AMINO_ACIDS},
                             {aa: float("nan") for aa in sr.AMINO_ACIDS}):
            modified = copy.deepcopy(model)
            modified["probabilities"]["1"] = distribution
            with self.assertRaises(ValueError):
                modified["model_sha256"] = sr._hash({k: v for k, v in modified.items() if k != "model_sha256"})
                sr.validate_profile(modified)

    def test_percentile_endpoints_and_midrank(self):
        rows = [record("a"), record("bb"), record("ccc", "C"),
                record("dddd", "C", "calibration"), record("eeeee", "A", "calibration")]
        model = fit(rows)
        for residue, expected in (("D", 0), ("C", 25), ("A", 75)):
            score = sr.score_numbered(model, {"1": residue, "2": "C"}, chain="H", domain_format="VH")
            self.assertEqual(score["percentile"], expected)

    def test_unknown_donor_explicitly_audited(self):
        for donor in ("unknown", "unknown:unit"):
            model = fit([record("a", donor=donor)])
            self.assertEqual(model["audit"]["unresolved_metadata"]["train"]["donor"], 1)
            self.assertEqual(model["audit"]["donor_support"]["train"], {"known_donors": 0, "proxy_donors": 1})
            self.assertEqual(model["evidence"], "research_only")

    def test_duplicate_conflicting_numbering_and_duplicate_ids(self):
        row = record("a")
        duplicate = dict(row, record_id="other", numbering={"1": "A", "3": "C", "128": "C"})
        with self.assertRaisesRegex(ValueError, "conflicting numbering"):
            fit([row, duplicate])
        with self.assertRaisesRegex(ValueError, "duplicate record_id"):
            fit([row, row])

    def test_exact_numbered_domain_required(self):
        row = record("a")
        for sequence in ("W" + row["sequence"], row["sequence"] + "W", "CAC"):
            with self.subTest(sequence=sequence), self.assertRaisesRegex(ValueError, "reconstructed numbered domain"):
                fit([dict(row, sequence=sequence)])
        other = record("bb", split="calibration")
        other["numbering"] = dict(row["numbering"])
        with self.assertRaisesRegex(ValueError, "reconstructed numbered domain"):
            fit([row, other])

    def test_imgt_order_independent_of_mapping_order(self):
        numbering = {"1": "A", "2": "C", "27": "D", "27A": "E", "33": "G", "33A": "F",
                     "61": "I", "61A": "H", "111": "K", "111A": "L", "111Z": "M",
                     "111AA": "N", "111AB": "P", "112": "A", "112A": "Y", "112Z": "W",
                     "112AA": "V", "112AB": "T", "128": "C"}
        row = record("a")
        row.update(numbering=numbering, sequence="ACDEFGHIKLMNPTVWYAC")
        model = fit([row])
        row["numbering"] = dict(reversed(list(numbering.items())))
        self.assertEqual(model, fit([row]))
        row["numbering"] = dict(sorted(numbering.items()))
        self.assertEqual(model, fit([row]))

    def test_mapping_hash_split_leakage_in_loaded_provenance(self):
        model = fit([record("a"), record("bb", split="calibration")])
        model["record_provenance"][1]["numbering_sha256"] = model["record_provenance"][0]["numbering_sha256"]
        model["model_sha256"] = sr._hash({k: v for k, v in model.items() if k != "model_sha256"})
        with self.assertRaisesRegex(ValueError, "provenance split leakage"):
            sr.validate_profile(model)

    def test_provenance_counts_revalidated(self):
        model = fit([record("a")])
        model["audit"]["effective_groups"]["train"]["records"] = 100
        model["model_sha256"] = sr._hash({k: v for k, v in model.items() if k != "model_sha256"})
        with self.assertRaisesRegex(ValueError, "provenance counts"):
            sr.validate_profile(model)

    def test_save_load_exclusive_and_fingerprint(self):
        model = fit([record("a")])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "model.json"
            sr.save_profile(model, path)
            self.assertEqual(sr.load_profile(path), model)
            with self.assertRaises(FileExistsError):
                sr.save_profile(model, path)
        audit = sr.audit_profile(model)
        self.assertFalse(audit["biologically_validated"])
        self.assertFalse(audit["provenance_verified"])
        audit["input_counts"]["train"] = 1000
        self.assertEqual(model["audit"]["input_counts"]["train"], 1)


if __name__ == "__main__":
    unittest.main()
