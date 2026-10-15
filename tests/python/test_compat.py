"""Compatibility contracts; synthetic inputs and alignment geometry only."""
import contextlib
import gzip
import importlib
import inspect
import io
from pathlib import Path
import tempfile
import unittest
import warnings

import anarci
import anarci_rs
from anarci import schemes
from anarci.germlines import all_germlines


class CompatibilityTests(unittest.TestCase):
    def test_signatures_and_defaults(self):
        self.assertEqual(list(inspect.signature(anarci.run_anarci).parameters), ["seq", "ncpu", "kwargs"])
        self.assertEqual(list(inspect.signature(anarci.number).parameters),
                         ["sequence", "scheme", "database", "allow", "allowed_species"])
        params = inspect.signature(anarci.anarci).parameters
        self.assertEqual(params["database"].default, "ALL")
        self.assertEqual(params["allowed_species"].default, ["human", "mouse"])
        self.assertEqual(inspect.signature(anarci_rs.anarci).parameters["database"].default, "pan")

    def test_imports_constants_and_germlines(self):
        scope = {"number": "caller binding"}
        exec("from anarci import *", scope)
        self.assertEqual(scope["number"], "caller binding")
        self.assertEqual(set(scope) - {"number", "__builtins__"}, {"anarci", "schemes"})
        self.assertTrue(callable(anarci.number))
        module = importlib.import_module("anarci.anarci")
        for name in ("run_hmmer", "parse_hmmer_output", "number_sequence_from_alignment",
                     "number_sequences_from_alignment", "get_identity", "run_germline_assignment",
                     "check_for_j", "_parse_hmmer_query", "_hmm_alignment_to_states", "validate_numbering"):
            self.assertTrue(callable(getattr(module, name)), name)
        self.assertEqual(schemes.alphabet[-1], " ")
        self.assertEqual(schemes.blosum62[("W", "W")], 11)
        self.assertIn("human", all_germlines["V"]["H"])
        self.assertEqual(anarci.chain_type_to_class["K"], "L")

    def test_silent_short_sequence_and_errors(self):
        with warnings.catch_warnings(record=True) as seen:
            self.assertEqual(anarci.number("ACDE"), (False, False))
            self.assertFalse(seen)
        with self.assertRaisesRegex(AssertionError, "Unknown amino acid"):
            anarci.number("A*")
        with self.assertRaisesRegex(AssertionError, "Unrecognised to unimplemented scheme"):
            anarci.number("ACDE", scheme="bad")
        with self.assertRaises(NotImplementedError):
            anarci.anarci([], hmmerpath="/custom/hmmscan")
        with self.assertRaisesRegex(AssertionError, "outfile"):
            anarci.anarci([], csv=True)

    def test_empty_batch_and_no_hit(self):
        self.assertEqual(anarci.anarci([]), ([], [], []))
        self.assertEqual(anarci.run_anarci([]), ([], [], [], []))
        result = anarci.anarci([("synthetic", "A" * 80)], allowed_species=None)
        self.assertEqual(result[0], [None])
        self.assertEqual(result[1], [None])
        self.assertEqual(result[2][0][0], ["id", "description", "evalue", "bitscore", "bias", "query_start", "query_end"])
        self.assertEqual(anarci.run_anarci("AAAA")[0], [["Input sequence", "AAAA"]])
        self.assertEqual(anarci.run_anarci([["synthetic", "AAAA"]])[0], [["synthetic", "AAAA"]])
        self.assertEqual(anarci.anarci([["synthetic", "AAAA"]])[0], [None])
        self.assertEqual(anarci.run_hmmer([["synthetic", "AAAA"]])[0][1], [])

    def test_empty_species_germline_legacy_semantics(self):
        from safe_reference import HEAVY
        with self.assertRaisesRegex(ValueError, r"max\(\) arg is an empty sequence"):
            anarci.anarci([("permitted-12e8-H", HEAVY)], assign_germline=True, allowed_species=[])
        result = anarci.anarci([("synthetic", "AAAA")], assign_germline=True, allowed_species=[])
        self.assertEqual(result[0], [None])
        result = anarci_rs.anarci([("permitted-12e8-H", HEAVY)], database="ALL",
                                  assign_germline=True, allowed_species=[])
        self.assertIn("germlines", result[1][0][0])

    def test_fasta_gzip_and_text_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            for suffix in (".fa", ".fa.gz"):
                path = Path(temp) / ("sequences" + suffix)
                opener = gzip.open if suffix.endswith(".gz") else open
                with opener(path, "wt") as stream:
                    stream.write(">synthetic full description\nACD\nEFG\n")
                self.assertEqual(anarci.read_fasta(str(path)), [("synthetic full description", "ACDEFG")])
                self.assertEqual(anarci.run_anarci(str(path))[0], [("synthetic full description", "ACDEFG")])
            text = io.StringIO()
            with contextlib.redirect_stdout(text):
                anarci.anarci([("synthetic", "AAAA")], output=True)
            self.assertEqual(text.getvalue(), "# synthetic\n//\n")
            output = Path(temp) / "numbering.txt"
            anarci.run_anarci([("synthetic", "AAAA")], output=True, outfile=output)
            self.assertEqual(output.read_text(), text.getvalue())

    def test_alignment_scheme_helpers(self):
        sequence = "A" * 128
        states = [((position, "m"), position - 1) for position in range(1, 129)]
        numbered = anarci.number_sequence_from_alignment(states, sequence)
        self.assertEqual(numbered[1:], (0, 127))
        nested_lists = [[[position, "m"], position - 1] for position in range(1, 129)]
        self.assertEqual(anarci.number_sequence_from_alignment(nested_lists, sequence), numbered)
        self.assertEqual(schemes.number_imgt(nested_lists, sequence), numbered)
        self.assertEqual("".join(aa for _, aa in numbered[0] if aa != "-"), sequence)
        for scheme, chain in (("imgt", None), ("aho", "H"), ("chothia", "H"),
                              ("kabat", "L"), ("martin", "H"), ("wolfguy", "L")):
            result = anarci.number_sequence_from_alignment(states, sequence, scheme, chain)
            self.assertIsInstance(result[0], list)
        self.assertEqual(anarci.get_identity(sequence, sequence), 1.)
        self.assertEqual(anarci.get_identity("-" * 128, "-" * 128), 0.)
        self.assertEqual(schemes.get_imgt_cdr(2, 12, 27, 39)[0], (27, " "))
        self.assertEqual(schemes.get_cdr3_annotations(13)[0], (105, " "))
        self.assertEqual(schemes.gap_missing([[((2, " "), "A")]]), [((1, " "), "-"), ((2, " "), "A")])

    def test_state_trust_boundary(self):
        for states in ([((0, "m"), 0)], [((129, "m"), 0)], [((1, "q"), 0)],
                       [((1, "m"), None)], [((1, "d"), 0)], [((1, "m"), 2)],
                       [((1, "m"), 0), ((2, "m"), 2)]):
            with self.subTest(states=states), self.assertRaises((AssertionError, ValueError)):
                anarci.number_sequence_from_alignment(states, "AA")
        with self.assertRaises(AssertionError):
            anarci.get_identity("A", "A")

    def test_hmmer_hsp_parser_without_biopython(self):
        from types import SimpleNamespace
        module = importlib.import_module("anarci.anarci")
        hsp = SimpleNamespace(hit_id="human_H", hit_description="synthetic", evalue=1e-20,
                              bitscore=100., bias=0., query_start=0, query_end=128,
                              hit_start=0, hit_end=128, order=0,
                              aln_annotation={"RF": "x" * 128, "PP": "*" * 128})
        result = module._parse_hmmer_query(SimpleNamespace(hsps=[hsp], seq_len=128))
        self.assertEqual(len(result[1][0]), 128)
        self.assertEqual(result[2][0]["chain_type"], "H")
        self.assertEqual(module._hmm_alignment_to_states(hsp, 1, 128), result[1][0])
        hsp.query_end = 129
        with self.assertRaises(AssertionError):
            module._hmm_alignment_to_states(hsp, 1, 128)

    def test_number_regions_private_helper(self):
        from anarci.schemes import _number_regions
        states = [((position, "m"), position - 1) for position in range(1, 129)]
        offsets = {0: 0}
        result = _number_regions("A" * 128, states, "X" * 128, "1" * 128, {"1": 0}, offsets, 1, [])
        self.assertEqual(len(result[0][0]), 128)
        self.assertEqual(result[1:], (0, 127))


if __name__ == "__main__":
    unittest.main()
