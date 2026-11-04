"""Numerical/reporting and pinned-source recipe contracts; no data downloads."""
import importlib.util
import json
import math
from pathlib import Path
import re
import unittest
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("species_evaluation", ROOT / "scripts/evaluate_species_reference.py")
evaluate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(evaluate)


class PipelineTests(unittest.TestCase):
    def test_auc_ties_empty_and_direction(self):
        self.assertEqual(evaluate.rank_auc([2, 3], [0, 1]), 1)
        self.assertEqual(evaluate.rank_auc([0, 1], [2, 3]), 0)
        self.assertEqual(evaluate.rank_auc([1, 1], [1, 1]), 0.5)
        self.assertIsNone(evaluate.rank_auc([], [1]))
        self.assertIsNone(evaluate.rank_auc([1], []))
        self.assertEqual(evaluate.rank_auc([1, 2], [0, 1]), 0.875)

    def test_fixed_panel_rejects_partial_score_comparisons(self):
        model = {"positions": ["1", "2"], "probabilities": {"1": {"A": 0.5}, "2": {"C": 0.25}}}
        self.assertAlmostEqual(evaluate.raw_score(model, {"1": "A", "2": "C"}), (math.log(.5) + math.log(.25)) / 2)
        self.assertIsNone(evaluate.raw_score(model, {"1": "A"}))
        self.assertIsNone(evaluate.raw_score(model, {"1": "A", "2": "X"}))

    def test_real_data_gate_rejects_synthetic_selected_or_unreviewed_sources(self):
        unit = {"url": "https://opig.stats.ox.ac.uk/example.csv.gz", "study": "study", "species": "human", "domain_format": "VH", "split": "train",
                "biological_origin": "natural_repertoire", "library_selection": "unselected",
                "origin_evidence": {"url": "https://example.org/paper", "review_note": "Numerical metadata contract fixture; not biological evidence"}}
        row = {**unit, "source_url": unit["url"]}
        evaluate.validate_cohort_origins([row], {"units": [unit]})
        for change in ({"biological_origin": "synthetic"}, {"library_selection": "selected"}, {"origin_evidence": {}}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                evaluate.validate_cohort_origins([row], {"units": [{**unit, **change}]})
        with self.assertRaises(ValueError):
            evaluate.validate_cohort_origins([{**row, "split": "test"}], {"units": [unit]})
        with self.assertRaises(ValueError):
            evaluate.validate_cohort_origins([{**row, "biological_origin": "synthetic"}], {"units": [unit]})

    def test_empty_distribution_does_not_invent_results(self):
        self.assertEqual(evaluate.distribution([]), {"n": 0, "minimum": None, "median": None, "maximum": None})
        self.assertEqual(evaluate.distribution([3, 1]), {"n": 2, "minimum": 1, "median": 2, "maximum": 3})

    def test_pinned_oas_recipe_has_real_source_integrity_and_separate_formats(self):
        manifest = json.loads((ROOT / ".agents/species/oas-sources.json").read_text())
        self.assertEqual(manifest["schema_version"], 1)
        self.assertEqual(len(manifest["units"]), 10)
        study_splits = {}
        species = set()
        total = 0
        for unit in manifest["units"]:
            parsed = urlsplit(unit["url"])
            self.assertEqual(parsed.scheme, "https")
            self.assertIn(parsed.hostname, {"opig.stats.ox.ac.uk", "oas-nano.opig.stats.ox.ac.uk"})
            self.assertTrue(parsed.path.endswith(".csv.gz"))
            self.assertRegex(unit["sha256"], r"^[0-9a-f]{64}$")
            self.assertGreater(unit["compressed_bytes"], 0)
            self.assertEqual(unit["biological_origin"], "natural_repertoire")
            self.assertEqual(unit["library_selection"], "unselected")
            self.assertTrue(unit["origin_evidence"]["url"])
            self.assertTrue(unit["origin_evidence"]["review_note"])
            self.assertNotIn(unit["study"], {"Chen_2021", "Rafique_2019", "Rafique_2023", "Jaffe_2022"})
            self.assertEqual(study_splits.setdefault(unit["study"], unit["split"]), unit["split"])
            self.assertIn(unit["split"], {"train", "calibration", "test"})
            species.add(unit["species"])
            total += unit["compressed_bytes"]
            if unit["domain_format"] == "VHH":
                self.assertTrue(unit["vhh_evidence"])
                self.assertEqual(parsed.hostname, "oas-nano.opig.stats.ox.ac.uk")
                self.assertEqual(unit["split"], "train")
        self.assertEqual(species, {"human", "mouse", "alpaca", "llama"})
        self.assertLess(total, 500_000_000)
        self.assertTrue(any("canine" in item for item in manifest["unsupported"]))
        self.assertEqual(manifest["selection"]["sample_size_per_unit"], 2000)


if __name__ == "__main__":
    unittest.main()
