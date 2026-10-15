"""Extended API boundary contracts; all input sequences are synthetic."""
import io
import math
from pathlib import Path
import tempfile
import unittest
import warnings

import anarci_rs


class ExtendedApiTests(unittest.TestCase):
    def test_extended_short_warning_and_none_species(self):
        with warnings.catch_warnings(record=True) as seen:
            warnings.simplefilter("always")
            self.assertEqual(anarci_rs.number("ACDE", allowed_species=None), (False, False))
            self.assertIn("min_length", str(seen[0].message))

    def test_extended_input_pairs_are_normalized(self):
        for records in ((("synthetic", "AAAA"),), [["synthetic", "AAAA"]],
                        iter([("synthetic", "AAAA")])):
            self.assertEqual(anarci_rs.run_anarci(records)[0], [("synthetic", "AAAA")])

    def test_native_scan_state_shapes(self):
        alignments = anarci_rs.run_hmmer([("synthetic", "A" * 80)], ncpu=2)
        self.assertEqual(len(alignments), 1)
        self.assertEqual(alignments[0][1:], ([], []))
        for sequence in ("A\x00A", "A-AA", "A☃A"):
            with self.subTest(sequence=sequence), self.assertRaises((AssertionError, ValueError)):
                anarci_rs.anarci([("synthetic", sequence)])
        with self.assertRaises((AssertionError, ValueError)):
            anarci_rs.run_hmmer([("synthetic", "AAAA")], bit_score_threshold=math.nan)

    def test_number_invalid_input_and_settings_raise(self):
        for sequence in ("A\x00A", "A" * 10000):
            with self.subTest(length=len(sequence)), self.assertRaises(AssertionError):
                anarci_rs.number(sequence)
        for threshold in (math.nan, math.inf, -math.inf):
            with self.subTest(threshold=threshold), self.assertRaises(AssertionError):
                anarci_rs.number("ACDE", bit_score_threshold=threshold)

    def test_version_and_native_features_exported(self):
        self.assertTrue(anarci_rs.__version__)
        for name in ("read_sequences", "translate_dna", "parse_structure", "annotate_protein",
                     "annotate_dna", "pair_structure", "format_annotations"):
            self.assertTrue(callable(getattr(anarci_rs, name)), name)

    def test_csv_output_helper(self):
        with tempfile.TemporaryDirectory() as temp:
            root = str(Path(temp) / "synthetic")
            numbered = [[([((1, " "), "A"), ((2, " "), "C")], 0, 1)]]
            details = [[dict(chain_type="K", species="human", evalue=1e-10, bitscore=100.)]]
            anarci_rs.csv_output([("synthetic,name", "AC")], numbered, details, root)
            lines = Path(root + "_KL.csv").read_text().splitlines()
            self.assertTrue(lines[0].endswith(",1,2"))
            self.assertTrue(lines[1].startswith("synthetic name,0,human,K"))
            self.assertTrue(lines[1].endswith(",A,C"))

    def test_write_and_validate_numbering(self):
        from anarci_rs.io import grouper
        stream = io.StringIO()
        anarci_rs.write_fasta([("synthetic", "A" * 81)], stream)
        self.assertEqual(stream.getvalue(), ">synthetic\n" + "A" * 80 + "\nA\n")
        self.assertEqual(list(grouper(2, range(3))), [[0, 1], [2]])
        numbered = ([((1, " "), "A"), ((2, " "), "C")], 0, 1)
        self.assertEqual(anarci_rs.validate_numbering(numbered, ("synthetic", "AC")), numbered)
        with self.assertRaises(AssertionError):
            anarci_rs.validate_numbering(numbered, ("synthetic", "AA"))


if __name__ == "__main__":
    unittest.main()
