"""Exact minimum-cardinality search over a caller-declared finite space.

Evidence and manifests are caller assertions, not authentication of callback
source or biological validation. Deadlines are cooperative, not preemptive.
Only the first accepted witness is returned; no best-score claim is made.
"""

import copy
import hashlib
import itertools
import json
import math
import time
from collections.abc import Sequence

__all__ = ["EvaluationError", "solve_minimum"]
_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


class EvaluationError(RuntimeError):
    """A callback failed or did not return its declared result contract."""

    def __init__(self, message, *, scope_manifest, scope_sha256, callback):
        super().__init__(message)
        self.scope_manifest = copy.deepcopy(scope_manifest)
        self.scope_sha256 = scope_sha256
        self.callback = callback


def _json_copy(value):
    # Reject lossy JSON conversions (non-string keys, tuples, NaN, etc.).
    def check(item):
        if item is None or type(item) in (str, bool, int):
            return
        if type(item) is float and math.isfinite(item):
            return
        if type(item) is list:
            for child in item:
                check(child)
            return
        if type(item) is dict and all(type(key) is str for key in item):
            for child in item.values():
                check(child)
            return
        raise ValueError("Expected finite JSON values with string object keys")
    check(value)
    return json.loads(json.dumps(value, sort_keys=True, allow_nan=False))


def _hash(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True, allow_nan=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def solve_minimum(parent, allowed, *, max_edits, screen, validators=(),
                  scope_metadata, max_nodes=200000, max_evaluations=100000,
                  max_validator_calls=10000, max_seconds=30.0, clock=None):
    """Return a JSON-compatible finite-scope proof, witness, or budget stop.

    Coordinates are arbitrary strings. Allowed alleles are canonical amino
    acids; wild type and duplicates are removed before enumeration. At most
    twelve positions may have non-wild-type alleles. Every position combination
    and every allele-product leaf costs one node (including both for k=0).
    A layer's evaluated count counts actual screens; rejected counts definitive
    screen/validator rejections. Validators run in order after screen acceptance.
    Counts at their exact caps do not prevent completion bookkeeping.

    Callbacks receive independent candidate copies and must return exactly
    {"accepted": bool, "evidence": dict}, with finite JSON evidence. Exceptions
    and malformed results raise EvaluationError, never scientific rejection.
    """
    if type(parent) is not dict or type(allowed) is not dict:
        raise ValueError("parent and allowed must be dictionaries")
    parent = copy.deepcopy(parent)
    allowed = copy.deepcopy(allowed)
    if any(type(key) is not str or type(aa) is not str or aa not in _AMINO_ACIDS
           for key, aa in parent.items()):
        raise ValueError("parent requires string coordinates and canonical amino acids")
    alleles = {}
    for key, values in allowed.items():
        if type(key) is not str or key not in parent:
            raise ValueError("allowed contains an unknown coordinate")
        if type(values) is not list or any(type(aa) is not str or aa not in _AMINO_ACIDS
                                          for aa in values):
            raise ValueError("allowed alleles must be lists of canonical amino acids")
        normalized = sorted(set(values) - {parent[key]})
        if normalized:
            alleles[key] = normalized
    if len(alleles) > 12:
        raise ValueError("At most 12 mutable positions are supported in total")
    if type(max_edits) is not int or not 0 <= max_edits <= 12:
        raise ValueError("max_edits must be an integer from 0 through 12")
    for name, value in (("max_nodes", max_nodes), ("max_evaluations", max_evaluations),
                        ("max_validator_calls", max_validator_calls)):
        if type(value) is not int or not 0 <= value <= 500000:
            raise ValueError(name + " must be an integer in 0..500000")
    if (type(max_seconds) not in (int, float) or not math.isfinite(max_seconds)
            or not 0 < max_seconds <= 300):
        raise ValueError("max_seconds must be finite and in (0, 300]")
    if not callable(screen) or (clock is not None and not callable(clock)):
        raise ValueError("screen and clock must be callable")
    if type(scope_metadata) is not dict:
        raise ValueError("scope_metadata must be a dictionary")
    metadata = _json_copy(scope_metadata)
    if not isinstance(validators, Sequence) or isinstance(validators, (str, bytes)):
        raise ValueError("validators must be a sequence of manifest/callback pairs")
    frozen_validators = []
    for entry in validators:
        if (not isinstance(entry, (tuple, list)) or len(entry) != 2
                or type(entry[0]) is not dict or not callable(entry[1])):
            raise ValueError("Each validator requires a manifest dictionary and callable")
        frozen_validators.append((_json_copy(entry[0]), entry[1]))
    scope = _json_copy({"parent": parent, "alleles": alleles, "max_edits": max_edits,
                        "metadata": metadata,
                        "validators": [entry[0] for entry in frozen_validators]})
    scope_hash = _hash(scope)
    budgets = {"max_nodes": max_nodes, "max_evaluations": max_evaluations,
               "max_validator_calls": max_validator_calls, "max_seconds": float(max_seconds)}
    run = {"scope_sha256": scope_hash, "budgets": budgets}
    counts = {"nodes": 0, "evaluations": 0, "validator_calls": 0}
    clock = time.monotonic if clock is None else clock
    start = clock()
    elapsed = 0.0
    excluded = None
    positions = sorted(alleles)
    layers = []
    for k in range(max_edits + 1):
        total = sum(math.prod(len(alleles[p]) for p in combo)
                    for combo in itertools.combinations(positions, k))
        layers.append({"k": k, "total_candidates": total, "evaluated": 0, "rejected": 0})

    class BudgetStop(Exception):
        pass

    def deadline():
        nonlocal elapsed
        elapsed = max(elapsed, float(clock() - start))
        if elapsed >= max_seconds:
            raise BudgetStop("max_seconds")

    def consume(counter, cap, label):
        deadline()
        if counts[counter] >= cap:
            raise BudgetStop(label)
        counts[counter] += 1

    def invoke(callback, candidate, label):
        try:
            result = callback(copy.deepcopy(candidate))
            if (type(result) is not dict or set(result) != {"accepted", "evidence"}
                    or type(result["accepted"]) is not bool
                    or type(result["evidence"]) is not dict):
                raise ValueError("Expected exactly accepted:boolean and evidence:object")
            result = _json_copy(result)
        except Exception as exc:
            raise EvaluationError("%s evaluation failed: %s" % (label, exc),
                                  scope_manifest=scope, scope_sha256=scope_hash,
                                  callback=label) from exc
        deadline()
        return result

    def finish(status, candidate=None, k=None, evidence=None, validation=None, reason=None):
        return {"status": status, "minimum_edits": k, "candidate": candidate,
                "screen_evidence": evidence, "validator_evidence": validation,
                "completed_excluded_through": excluded, "layers": layers,
                "search": dict(counts, elapsed_seconds=elapsed),
                "scope_manifest": scope, "scope_sha256": scope_hash,
                "run_manifest": run, "run_sha256": _hash(run),
                "proof": status != "budget_exhausted",
                "proof_scope": "finite_declared_space", "biologically_validated": False,
                "budget_exhausted_reason": reason,
                "deadline_overrun_seconds": max(0.0, elapsed - max_seconds)}

    try:
        for layer in layers:
            k = layer["k"]
            for combo in itertools.combinations(positions, k):
                consume("nodes", max_nodes, "max_nodes")
                for replacements in itertools.product(*(alleles[p] for p in combo)):
                    consume("nodes", max_nodes, "max_nodes")
                    candidate = dict(parent)
                    candidate.update(zip(combo, replacements))
                    consume("evaluations", max_evaluations, "max_evaluations")
                    layer["evaluated"] += 1
                    result = invoke(screen, candidate, "screen")
                    accepted = result["accepted"]
                    validation = []
                    if accepted:
                        for index, (manifest, callback) in enumerate(frozen_validators):
                            consume("validator_calls", max_validator_calls, "max_validator_calls")
                            checked = invoke(callback, candidate, "validator[%d]" % index)
                            validation.append({"manifest": copy.deepcopy(manifest),
                                               "evidence": checked["evidence"]})
                            if not checked["accepted"]:
                                accepted = False
                                break
                    if accepted:
                        return finish("parent_meets_target" if k == 0 else "proven_minimum",
                                      candidate, k, result["evidence"], validation)
                    layer["rejected"] += 1
            excluded = k
        return finish("proven_infeasible_within_scope")
    except BudgetStop as exc:
        return finish("budget_exhausted", reason=str(exc))
