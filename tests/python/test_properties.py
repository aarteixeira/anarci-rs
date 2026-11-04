"""Bounded generated contracts for the installed compatibility and extended APIs.

No historical fixtures are opened. Native integration uses only the existing
safe_reference public 12e8 sequences and explicitly synthetic derivatives.
"""
import inspect
import math
import unittest
import warnings
from unittest import mock

import anarci
import anarci_rs
from anarci_rs import _api
from safe_reference import HEAVY, LIGHT


class FacadePropertyTests(unittest.TestCase):
    def test_separate_defaults_are_forwarded_without_cross_facade_mutation(self):
        # Capture the actual call, rather than checking signatures alone.
        for facade, database, species in ((anarci, "ALL", ["human", "mouse"]),
                                           (anarci_rs, "pan", None)):
            for _ in range(2):
                with self.subTest(facade=facade.__name__), mock.patch.object(
                        _api._native, "anarci", return_value=([], [], [])) as native:
                    self.assertEqual(facade.anarci([]), ([], [], []))
                    args, kwargs = native.call_args
                    self.assertEqual(args, ([],))
                    self.assertEqual(kwargs["database"], database)
                    self.assertEqual(kwargs["allowed_species"], species)
                    self.assertEqual(kwargs["scheme"], "imgt")
                    self.assertFalse(kwargs["assign_germline"])
                    self.assertEqual(kwargs["bit_score_threshold"], 80)
        self.assertNotIn("min_length", inspect.signature(anarci.number).parameters)
        self.assertEqual(inspect.signature(anarci_rs.number).parameters["min_length"].default, 70)

    def test_facade_option_errors_precede_native_calls_even_for_empty_batches(self):
        for facade in (anarci, anarci_rs):
            for kwargs, error in (({"scheme": "not-a-scheme"}, AssertionError),
                                  ({"csv": True}, AssertionError),
                                  ({"hmmerpath": "/not/an/authorized/executable"}, NotImplementedError)):
                with self.subTest(facade=facade.__name__, kwargs=kwargs), mock.patch.object(
                        _api._native, "anarci") as native:
                    with self.assertRaises(error):
                        facade.anarci([], **kwargs)
                    native.assert_not_called()

    def test_short_sequence_warning_contracts_stay_separate(self):
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        for length in (0, 1, 2, 19, 69):
            sequence = (alphabet * 4)[:length]
            with self.subTest(length=length), warnings.catch_warnings(record=True) as seen:
                warnings.simplefilter("always")
                self.assertEqual(anarci.number(sequence), (False, False))
                self.assertEqual(seen, [])
                self.assertEqual(anarci_rs.number(sequence), (False, False))
                self.assertEqual(len(seen), 1)
                self.assertIn("min_length", str(seen[0].message))

    def test_invalid_states_are_python_errors_for_both_facades(self):
        for facade in (anarci, anarci_rs):
            for position in (0, 129, 255):
                with self.subTest(facade=facade.__name__, position=position):
                    with self.assertRaises((AssertionError, ValueError, OverflowError)):
                        facade.number_sequence_from_alignment([((position, "m"), 0)], "A")
            for states in ([((23, "m"), None)], [((23, "d"), 0)],
                           [((23, "m"), 1)], [((23, "?"), 0)],
                           [((23, "m"), 0), ((24, "m"), 0)]):
                with self.subTest(facade=facade.__name__, states=states):
                    with self.assertRaises((AssertionError, ValueError, OverflowError)):
                        facade.number_sequence_from_alignment(states, "A")

    def test_generated_alignment_geometry_preserves_subsequence_and_coordinates(self):
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        sequence = "GG" + (alphabet * 7)[:128] + "SS"
        for facade in (anarci, anarci_rs):
            for offset in (0, 1, 2):
                states = [((position, "m"), position - 1 + offset) for position in range(1, 129)]
                for scheme, alias in (("imgt", "IMGT"), ("kabat", "KABAT"),
                                      ("chothia", "CHOTHIA"), ("martin", "MARTIN"),
                                      ("aho", "AHO"), ("wolfguy", "WOLFGUY")):
                    with self.subTest(facade=facade.__name__, offset=offset, scheme=scheme):
                        result = facade.number_sequence_from_alignment(states, sequence, scheme, "H")
                        self.assertEqual(result, facade.number_sequence_from_alignment(states, sequence, alias, "H"))
                        self.assertEqual(result[1:], (offset, offset + 127))
                        self.assertEqual("".join(aa for _, aa in result[0] if aa != "-"),
                                         sequence[offset:offset + 128])

    def test_scan_input_and_nonfinite_settings_fail_instead_of_no_hit(self):
        for facade in (anarci, anarci_rs):
            for name, sequence in (("nul\0name", "AAA"), ("synthetic", "AA-AA"),
                                   ("synthetic", "A\0A"), ("synthetic", "A☃A")):
                with self.subTest(facade=facade.__name__, name=name, sequence=sequence):
                    with self.assertRaises((AssertionError, ValueError)):
                        facade.anarci([("valid", "AAAA"), (name, sequence)], database="ALL")
            for threshold in (math.nan, math.inf, -math.inf):
                with self.subTest(facade=facade.__name__, threshold=threshold):
                    with self.assertRaises((AssertionError, ValueError)):
                        facade.anarci([], database="ALL", bit_score_threshold=threshold)

    def test_installed_native_batch_matches_singletons_and_restores_query_metadata(self):
        # Public safe helper only; suffixes are synthetic and outside the V domain.
        records = [("heavy-first", HEAVY), ("no-hit", "A" * 80), ("light", LIGHT),
                   ("heavy-copy", HEAVY), ("heavy-suffix", HEAVY + "ACDE"), ("light-copy", LIGHT)]
        for facade in (anarci, anarci_rs):
            singles = [facade.anarci([record], database="ALL", allowed_species=None,
                                      ncpu=1) for record in records]
            for workers in (1, 2):
                result = facade.anarci(records, database="ALL", allowed_species=None, ncpu=workers)
                self.assertEqual(len(result), 3)
                for column in range(3):
                    self.assertEqual(result[column], [single[column][0] for single in singles])
                self.assertIsNone(result[0][1])
                self.assertIsNone(result[1][1])
                for index in (0, 2, 3, 4, 5):
                    self.assertTrue(result[0][index], "safe public input unexpectedly produced no domain")
                    for detail in result[1][index]:
                        self.assertEqual(detail["query_name"], records[index][0])
                        self.assertEqual(detail["scheme"], "imgt")
                result[1][0][0]["query_name"] = "mutated by caller"
                self.assertEqual(result[1][3][0]["query_name"], "heavy-copy")


if __name__ == "__main__":
    unittest.main()
