"""Synthetic finite-space contracts; no native extension or biological claims."""
import importlib.util
import itertools
import json
from pathlib import Path
import random
import unittest

_SOURCE = Path(__file__).resolve().parents[2] / "python" / "anarci_rs" / "minimum_edit.py"
if _SOURCE.is_file():
    _SPEC = importlib.util.spec_from_file_location("standalone_minimum_edit", _SOURCE)
    _MODULE = importlib.util.module_from_spec(_SPEC)
    _SPEC.loader.exec_module(_MODULE)
    solve_minimum = _MODULE.solve_minimum
    EvaluationError = _MODULE.EvaluationError
else:
    from anarci_rs.minimum_edit import EvaluationError, solve_minimum


def answer(accepted, **evidence):
    return {"accepted": accepted, "evidence": evidence}


class MinimumEditTests(unittest.TestCase):
    def solve(self, parent=None, allowed=None, **kwargs):
        options = dict(max_edits=2, screen=lambda c: answer(False), scope_metadata={})
        options.update(kwargs)
        return solve_minimum({"H:42": "A"} if parent is None else parent,
                             {"H:42": ["C"]} if allowed is None else allowed, **options)

    def test_empty_parent_zero_witness(self):
        result = self.solve({}, {}, max_edits=0, screen=lambda c: answer(True))
        self.assertEqual(result["status"], "parent_meets_target")
        self.assertEqual(result["candidate"], {})
        self.assertEqual(result["minimum_edits"], 0)
        self.assertIsNone(result["completed_excluded_through"])
        self.assertEqual(result["search"]["nodes"], 2)
        self.assertTrue(result["proof"])
        json.dumps(result, allow_nan=False)

    def test_max_zero_and_empty_alleles(self):
        for allowed in ({}, {"H:42": ["A", "A"]}):
            result = self.solve(allowed=allowed, max_edits=0)
            self.assertEqual(result["status"], "proven_infeasible_within_scope")
            self.assertEqual(result["completed_excluded_through"], 0)
            self.assertEqual(result["search"]["evaluations"], 1)

    def test_independent_cartesian_oracle(self):
        rng = random.Random(2718)
        parent = {"H:42": "A", "L:7": "C", "arbitrary": "D"}
        allowed = {"H:42": ["C", "D"], "L:7": ["A"], "arbitrary": ["E", "F"]}
        keys = list(parent)
        universe = [dict(zip(keys, aas)) for aas in itertools.product(
            *[[parent[p]] + allowed[p] for p in keys])]
        for _ in range(40):
            accepted = {tuple(c.items()) for c in universe if rng.random() < 0.15}
            cap = rng.randrange(4)
            distances = [sum(c[p] != parent[p] for p in keys) for c in universe
                         if tuple(c.items()) in accepted]
            feasible = [d for d in distances if d <= cap]
            result = self.solve(parent, allowed, max_edits=cap,
                                screen=lambda c: answer(tuple(c.items()) in accepted))
            self.assertEqual(result["minimum_edits"], min(feasible) if feasible else None)
            self.assertEqual(result["status"], ("parent_meets_target" if min(feasible) == 0
                             else "proven_minimum") if feasible
                             else "proven_infeasible_within_scope")
            if feasible:
                k = min(feasible)
                self.assertEqual(result["completed_excluded_through"], k - 1 if k else None)
                for layer in result["layers"][:k]:
                    self.assertEqual(layer["rejected"], layer["total_candidates"])

    def test_compensating_guard_and_motif_repair_paired_total(self):
        # Single edits fail a joint final guard; one edit lowers a toy score.
        parent = {"H:42": "A", "L:42": "C"}
        allowed = {"H:42": ["D"], "L:42": ["E"]}
        visited = []
        def screen(c):
            visited.append(dict(c))
            score = (2 if c["H:42"] == "D" else 0) - (1 if c["L:42"] == "E" else 0)
            return answer(score >= 1)
        def guard(c):
            return answer(c["H:42"] == "D" and c["L:42"] == "E", motif="repaired")
        result = self.solve(parent, allowed, screen=screen, validators=[({"id": "joint"}, guard)])
        self.assertEqual(result["minimum_edits"], 2)
        self.assertEqual(len(visited), 4)
        self.assertEqual(result["validator_evidence"][0]["evidence"], {"motif": "repaired"})
        self.assertFalse(self.solve(parent, allowed, max_edits=1, screen=screen,
                                   validators=[({}, guard)])["candidate"])

    def test_exact_caps_complete_without_phantom_nodes(self):
        result = self.solve(max_edits=1, max_nodes=4, max_evaluations=2)
        self.assertEqual(result["status"], "proven_infeasible_within_scope")
        self.assertEqual(result["search"]["nodes"], 4)
        for cap in range(4):
            result = self.solve(max_edits=1, max_nodes=cap)
            self.assertEqual(result["status"], "budget_exhausted")
            self.assertLessEqual(result["search"]["nodes"], cap)
            self.assertFalse(result["proof"])
        for cap in (0, 1):
            result = self.solve(max_evaluations=cap)
            self.assertEqual(result["search"]["evaluations"], cap)
            self.assertEqual(result["status"], "budget_exhausted")
        result = self.solve(max_edits=0, max_validator_calls=1, screen=lambda c: answer(True),
                            validators=[({}, lambda c: answer(False))])
        self.assertEqual(result["status"], "proven_infeasible_within_scope")
        result = self.solve(max_validator_calls=1, screen=lambda c: answer(True),
                            validators=[({}, lambda c: answer(True)), ({}, lambda c: answer(True))])
        self.assertEqual(result["status"], "budget_exhausted")
        self.assertEqual(result["search"]["validator_calls"], 1)
        self.assertIsNone(result["candidate"])

    def test_cooperative_deadline_overrun(self):
        now = [0.0]
        def slow(c):
            now[0] += 2
            return answer(True)
        result = self.solve(clock=lambda: now[0], max_seconds=1, screen=slow)
        self.assertEqual(result["status"], "budget_exhausted")
        self.assertEqual(result["budget_exhausted_reason"], "max_seconds")
        self.assertEqual(result["deadline_overrun_seconds"], 1)
        self.assertEqual(result["search"]["evaluations"], 1)
        self.assertIsNone(result["minimum_edits"])
        self.assertEqual(result["layers"][0]["rejected"], 0)

    def test_bad_callbacks_raise_typed_error(self):
        bad = [None, {}, {"accepted": 1, "evidence": {}},
               {"accepted": True, "evidence": []},
               {"accepted": True, "evidence": {"x": float("nan")}},
               {"accepted": False, "evidence": {}, "unknown": 1}]
        def fail(c):
            raise RuntimeError("native failed")
        for cb in [lambda c, value=v: value for v in bad] + [fail]:
            for validator in (False, True):
                with self.assertRaises(EvaluationError) as raised:
                    self.solve(screen=(lambda c: answer(True)) if validator else cb,
                               validators=[({}, cb)] if validator else ())
                self.assertEqual(len(raised.exception.scope_sha256), 64)
                self.assertEqual(raised.exception.scope_manifest["parent"], {"H:42": "A"})

    def test_hashes_normalization_and_budget_separation(self):
        first = self.solve(allowed={"H:42": ["A", "D", "C", "C"]})
        second = self.solve(allowed={"H:42": ["C", "D"]}, max_nodes=999)
        self.assertEqual(first["scope_sha256"], second["scope_sha256"])
        self.assertNotEqual(first["run_sha256"], second["run_sha256"])
        third = self.solve(allowed={"H:42": ["D", "C"]})
        self.assertEqual(first["run_sha256"], third["run_sha256"])
        self.assertNotEqual(first["scope_sha256"], self.solve(scope_metadata={"source": "x"})["scope_sha256"])

    def test_inputs_and_witness_are_frozen(self):
        parent = {"H:42": "A"}
        allowed = {"H:42": ["C"]}
        metadata = {"nested": [1]}
        manifest = {"nested": [2]}
        evidence = {"nested": [3]}
        def screen(c):
            c["H:42"] = "Y"
            parent["H:42"] = "Y"
            allowed["H:42"].append("Y")
            metadata["nested"].append(9)
            manifest["nested"].append(9)
            return {"accepted": True, "evidence": evidence}
        def validator(c):
            self.assertEqual(c, {"H:42": "A"})
            c.clear()
            evidence["nested"].append(99)
            return answer(True)
        result = self.solve(parent, allowed, screen=screen, scope_metadata=metadata,
                            validators=[(manifest, validator)])
        self.assertEqual(result["candidate"], {"H:42": "A"})
        self.assertEqual(result["scope_manifest"]["metadata"], {"nested": [1]})
        self.assertEqual(result["scope_manifest"]["validators"], [{"nested": [2]}])
        self.assertEqual(result["screen_evidence"], {"nested": [3]})
        # Solver itself does not mutate any inputs.
        parent, allowed = {"x": "A"}, {"x": ["C", "A", "C"]}
        self.solve(parent, allowed)
        self.assertEqual(parent, {"x": "A"})
        self.assertEqual(allowed, {"x": ["C", "A", "C"]})

    def test_invalid_inputs(self):
        for options in [dict(max_edits=True), dict(max_edits=-1), dict(max_edits=13),
                        dict(max_nodes=True), dict(max_evaluations=-1),
                        dict(max_validator_calls=1.5), dict(max_seconds=True),
                        dict(max_seconds=0), dict(max_seconds=301),
                        dict(max_seconds=float("inf")), dict(max_seconds=float("nan")),
                        dict(scope_metadata={"x": float("nan")}), dict(scope_metadata={1: "x"}),
                        dict(screen=None), dict(validators=[({}, None)])]:
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.solve(**options)
        for parent, allowed in [({"x": True}, {}), ({"x": "X"}, {}),
                                ({1: "A"}, {}), ({"x": "A"}, {"y": ["C"]}),
                                ({"x": "A"}, {"x": [False]}),
                                ({"x": "A"}, {"x": "C"}),
                                ({str(i): "A" for i in range(13)},
                                 {str(i): ["C"] for i in range(13)})]:
            with self.assertRaises(ValueError):
                self.solve(parent, allowed)


if __name__ == "__main__":
    unittest.main()
