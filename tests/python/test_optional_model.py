"""CPU research model smoke checks; only short synthetic strings are used."""
import math
import unittest

from anarci_rs import experimental as exp


class ModelAvailabilityTests(unittest.TestCase):
    @unittest.skipIf(exp.MODEL_AVAILABLE, "model feature enabled")
    def test_default_build_explains_optional_feature(self):
        for call in (lambda: exp.Model(), lambda: exp.tokenize("AC"),
                     lambda: exp.overfit([("AC", ["out", "out"])])):
            with self.assertRaisesRegex(RuntimeError, "learned-model"):
                call()


@unittest.skipUnless(exp.MODEL_AVAILABLE, "requires learned-model feature")
class OptionalModelTests(unittest.TestCase):
    def test_forward_and_logits_shapes(self):
        model = exp.Model.new(dim=4, layers=1, ff=8)
        tokens = exp.tokenize("AcX")
        self.assertEqual(tokens, [1, 2, 21])
        logits = model.logits("AC")
        self.assertEqual((len(logits), len(logits[0])), (2, exp.NUM_LABELS))
        self.assertTrue(all(math.isfinite(score) for row in logits for score in row))
        batch = model.forward([tokens, tokens])
        self.assertEqual((len(batch), len(batch[0]), len(batch[0][0])), (2, 3, exp.NUM_LABELS))

    def test_one_training_step(self):
        model, accuracy = exp.overfit([("AC", [("match", 1), ("match", 2)])],
                                     dim=4, layers=1, ff=8, steps=1, lr=0.001)
        self.assertIsInstance(model, exp.Model)
        self.assertTrue(math.isfinite(accuracy))
        self.assertTrue(0.0 <= accuracy <= 1.0)
        self.assertEqual(len(model.logits("AC")), 2)

    def test_invalid_model_inputs(self):
        with self.assertRaises(ValueError):
            exp.Model(dim=0)
        model = exp.Model(dim=4, layers=0, ff=8)
        for tokens in ([], [[]], [[22]], [[1], [1, 2]]):
            with self.assertRaises(ValueError):
                model.forward(tokens)
        for sequence in ("", "é", "A" * 1025):
            with self.assertRaises(ValueError):
                model.logits(sequence)
        for examples in ([], [("AC", ["out"])]):
            with self.assertRaises(ValueError):
                exp.overfit(examples, dim=4, layers=0, ff=8)
        with self.assertRaises(ValueError):
            exp.overfit([("AC", ["out", "out"])], lr=float("nan"))


if __name__ == "__main__":
    unittest.main()
