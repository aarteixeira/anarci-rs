"""Research API tests using short, synthetic inputs only."""
import unittest

from anarci_rs import experimental as exp


class ExperimentalTests(unittest.TestCase):
    def test_label_vocabulary_round_trip(self):
        for index in range(exp.NUM_LABELS):
            self.assertEqual(exp.label_to_index(exp.index_to_label(index)), index)
        for invalid in (-1, 261, True, 1.5):
            with self.assertRaises(ValueError):
                exp.index_to_label(invalid)
        for invalid in (("match", 0), ("insert", 129), "unknown", ("match", True)):
            with self.assertRaises(ValueError):
                exp.label_to_index(invalid)

    def test_codec_restores_deletes_and_insertions(self):
        states = [(1, "M", 0), (2, "D", None), (3, "I", 1), (3, "M", 2)]
        encoded = exp.encode_domain(states)
        self.assertEqual(encoded, [(0, ("match", 1)), (1, ("insert", 3)), (2, ("match", 3))])
        self.assertEqual(exp.decode_domain(encoded), states)
        self.assertEqual(exp.strip_boundary_deletes([(1, "D", None), (2, "M", 0), (3, "D", None)]), [(2, "M", 0)])
        self.assertEqual(exp.decode_sequence([("match", 1), ("match", 2), "linker", ("match", 1)]),
                         [[(1, "M", 0), (2, "M", 1)], [(1, "M", 3)]])

    def test_bad_states_raise_python_exceptions(self):
        for states in (
            [(0, "M", 0)], [(129, "M", 0)], [(1, "?", 0)],
            [(1, "D", 0)], [(1, "M", None)],
            [(1, "M", 0), (1, "M", 1)],
            [(1, "M", 1), (2, "M", 0)],
        ):
            with self.subTest(states=states), self.assertRaises((ValueError, OverflowError)):
                exp.encode_domain(states)
        with self.assertRaises(ValueError):
            exp.DomainRecord("AC", [(1, "M", 2)])
        with self.assertRaises(ValueError):
            exp.decode_domain([(0, ("match", 2)), (1, ("match", 1))])

    def test_constrained_decode_and_score_validation(self):
        logits = [[-20.0] * exp.NUM_LABELS for _ in range(3)]
        for row, label in zip(logits, [("match", 1), ("insert", 2), ("match", 2)]):
            row[exp.label_to_index(label)] = 20.0
        self.assertEqual(exp.viterbi_decode(logits), [("match", 1), ("insert", 2), ("match", 2)])
        self.assertEqual(exp.viterbi_decode([]), [])
        for bad in ([[0.0]], [[0.0] * 262], [[float("nan")] * 261], [[float("inf")] * 261]):
            with self.assertRaises(ValueError):
                exp.viterbi_decode(bad)

    def test_augmentation_preserves_labels_and_boundaries(self):
        segment = exp.LabeledSeq.from_domain("ACD", [(1, "M", 0), (2, "M", 1), (3, "M", 2)])
        self.assertEqual(segment.restrict_positions(2, 3).sequence, "CD")
        self.assertEqual(segment.restrict_window(1, 99).sequence, "CD")
        self.assertEqual(segment.restrict_window(9, 20).sequence, "")
        self.assertEqual(segment.domains(), [[(1, "M", 0), (2, "M", 1), (3, "M", 2)]])
        self.assertEqual(exp.LabeledSeq.junk("AA").labels, ("out", "out"))
        self.assertEqual(exp.LabeledSeq.region("GG", "hinge").labels, ("hinge", "hinge"))
        self.assertEqual(exp.gly_ser_linker(2).sequence, "GGGGSGGGGS")
        self.assertEqual(exp.gly_ser_linker(0).sequence, "")
        self.assertEqual(len(exp.concat([segment, segment]).domains()), 2)
        joined = exp.concat_domains([segment, segment])
        self.assertEqual(joined.sequence, "ACDGACD")
        self.assertEqual(len(joined.domains()), 2)
        with self.assertRaises(ValueError):
            exp.LabeledSeq("AC", ("out",))
        with self.assertRaises(ValueError):
            exp.LabeledSeq.junk("é")
        with self.assertRaises(ValueError):
            exp.LabeledSeq.region("A", ("match", 1))
        with self.assertRaises(ValueError):
            segment.restrict_positions(3, 2)
        with self.assertRaises((ValueError, OverflowError)):
            exp.gly_ser_linker(2 ** 63)

    def test_dataset_helpers(self):
        a = exp.DomainRecord("ACDE", [(105, "M", 0), (106, "M", 1)])
        b = exp.DomainRecord("ACDF", [(105, "M", 0), (106, "M", 1)])
        c = exp.DomainRecord("WWWW", [(105, "M", 0), (106, "M", 1)])
        self.assertEqual(exp.cdr3_of(a), "AC")
        self.assertEqual(exp.overlap_identity("AC", "TAC"), 1.0)
        self.assertEqual(exp.overlap_identity("", "AC"), 0.0)
        self.assertTrue(exp.near_duplicate(a, b))
        self.assertFalse(exp.near_duplicate(a, c))
        clusters = exp.cluster_records([a, b, c])
        self.assertEqual(clusters[0], clusters[1])
        self.assertNotEqual(clusters[0], clusters[2])
        splits = exp.splits_of([a, b, c])
        self.assertEqual(splits[0], splits[1])
        self.assertIn(exp.split_for_rep("ACDE"), ("train", "val", "test"))
        examples = exp.build_dataset([a, b], [(105, 106)])
        self.assertTrue({"single", "partial", "multi_v", "junk"}.issubset({e["kind"] for e in examples}))
        for example in examples:
            self.assertIsInstance(example["data"], exp.LabeledSeq)
            self.assertEqual(len(example["data"].sequence), len(example["data"].labels))
        with self.assertRaises(ValueError):
            exp.cluster_records([a], full_thr=float("nan"))
        with self.assertRaises(ValueError):
            exp.cluster_records([a], cdr3_thr=1.1)
        with self.assertRaises(ValueError):
            exp.build_dataset([a], [(0, 128)])


if __name__ == "__main__":
    unittest.main()
