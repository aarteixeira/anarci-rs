"""Synthetic mapping contracts only; these coordinates are not biological evidence."""
import copy
import hashlib
import json
import unittest
from unittest.mock import patch

from anarci_rs import autonomy_structure as module


def chain(name, sequence, coordinates):
    return {"id": name, "sequence": sequence,
            "residues": [{"aa": aa, "ca": ca} for aa, ca in zip(sequence, coordinates)]}


def numbered(sequence, domain_format):
    start, end = (1, 3) if domain_format == "VH" else (1, 2)
    return {"chain": "H" if domain_format == "VH" else "K", "domain_count": 1,
            "start": start, "end": end,
            "numbering": {str(i + 40): sequence[i] for i in range(start, end + 1)},
            "sequence_indices": {str(i + 40): i for i in range(start, end + 1)}}


class InterfaceEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.chains = [chain("12E8_H", "WACDF", [(99, 0, 0), (0, 0, 0), (20, 0, 0), (40, 0, 0), (99, 0, 0)]),
                       chain("12E8_L", "QEFG", [(40, 0, 0), (0, 8, 0), (20, 9, 0), (40, 0, 0)]),
                       chain("12E8_A", "Y", [(20, 10, 0)])]
        self.provenance = {"source_id": "synthetic-contract", "structure_kind": "predicted",
                           "assembly_context": "unknown", "reviewed": True}
        self.parse = patch.object(module, "parse_structure", return_value=self.chains).start()
        self.number = patch.object(module, "number_domain", side_effect=numbered).start()
        self.addCleanup(patch.stopall)

    def build(self, **kwargs):
        args = dict(vh_chain_id="12E8_H", vl_chain_id="12E8_L", provenance=self.provenance)
        args.update(kwargs)
        return module.build_interface_evidence("GACDT", "synthetic structure\n", **args)

    def test_distances_mapping_scope_hash_and_json(self):
        result = self.build(antigen_chain_ids=["12E8_A"])
        self.assertEqual(result["interface_positions"], ["41"])
        self.assertEqual(result["antigen_protected_positions"], ["42"])
        self.assertEqual(result["position_evidence"][0], {
            "position": "41", "sequence_index": 1, "chain_residue_index": 1,
            "min_vl_ca_distance": 8.0, "min_antigen_ca_distance": 500 ** 0.5})
        self.assertGreater(result["position_evidence"][2]["min_vl_ca_distance"], 8)
        self.assertEqual(result["query_sequence_sha256"], hashlib.sha256(b"GACDT").hexdigest())
        self.assertEqual(result["structure_sha256"], hashlib.sha256(b"synthetic structure\n").hexdigest())
        self.assertEqual(result["structure_assessment"], "parent_CA_proximity_only")
        self.assertFalse(result["experimental_validation"])
        self.assertEqual(result["interface_scope"], "frozen_parent")
        self.assertNotIn("exposed", result)
        self.assertNotIn("binding_preserved", result)
        json.dumps(result, allow_nan=False)
        self.number.assert_any_call("GACDT", domain_format="VH")
        self.number.assert_any_call("QEFG", domain_format="VL")

    def test_consumed_evidence_cannot_drop_antigen_locks_or_required_schema(self):
        evidence = self.build(antigen_chain_ids=["12E8_A"])
        module.validate_interface_evidence(evidence, "GACDT", evidence["domain"])
        for key in ("chain_roles", "metrics_definition", "position_evidence", "provenance"):
            changed = copy.deepcopy(evidence); changed.pop(key)
            with self.subTest(key=key), self.assertRaises(ValueError):
                module.validate_interface_evidence(changed, "GACDT", evidence["domain"])
        changed = copy.deepcopy(evidence); changed["antigen_protected_positions"] = []
        with self.assertRaisesRegex(ValueError, "antigen_protected_positions"):
            module.validate_interface_evidence(changed, "GACDT", evidence["domain"])
        changed = copy.deepcopy(evidence); changed["interface_positions"] = ["41", "42"]
        with self.assertRaisesRegex(ValueError, "interface_positions"):
            module.validate_interface_evidence(changed, "GACDT", evidence["domain"])

    def test_absent_antigen_is_unknown(self):
        result = self.build()
        self.assertIsNone(result["position_evidence"][0]["min_antigen_ca_distance"])
        self.assertIn("antigen proximity unknown", result["unknowns"][0])
        self.assertEqual(result["antigen_protected_positions"], [])

    def test_exact_ids_and_distinct_roles(self):
        for kwargs in ({"vh_chain_id": "H"}, {"vl_chain_id": "missing"},
                       {"vl_chain_id": "12E8_H"}, {"antigen_chain_ids": ["missing"]},
                       {"antigen_chain_ids": ["12E8_H"]}, {"antigen_chain_ids": ["12E8_L"]},
                       {"antigen_chain_ids": ["12E8_A", "12E8_A"]}, {"antigen_chain_ids": "12E8_A"}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                self.build(**kwargs)

    def test_duplicate_native_id(self):
        self.chains.append(copy.deepcopy(self.chains[0]))
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            self.build()

    def test_missing_and_ambiguous_domain_substring(self):
        for observed in ("WADF", "WACDACDF"):
            self.chains[0] = chain("12E8_H", observed, [(0, 0, 0)] * len(observed))
            with self.subTest(observed=observed), self.assertRaisesRegex(ValueError, "unique exact"):
                self.build()

    def test_missing_nonfinite_or_malformed_coordinates_fail(self):
        for chain_index, residue_index in ((0, 1), (0, 3), (1, 1), (1, 2), (2, 0)):
            for ca in (None, [None, 0, 0], [float("nan"), 0, 0], [float("inf"), 0, 0], [0, 0], [True, 0, 0]):
                old = self.chains[chain_index]["residues"][residue_index]["ca"]
                self.chains[chain_index]["residues"][residue_index]["ca"] = ca
                with self.subTest(chain=chain_index, ca=ca), self.assertRaisesRegex(ValueError, "CA coordinate"):
                    self.build(antigen_chain_ids=["12E8_A"])
                self.chains[chain_index]["residues"][residue_index]["ca"] = old

    def test_constant_region_coordinates_not_used(self):
        self.chains[0]["residues"][0]["ca"] = None
        self.chains[1]["residues"][0]["ca"] = None
        self.assertEqual(self.build()["interface_positions"], ["41"])

    def test_no_interface(self):
        with self.assertRaisesRegex(ValueError, "No VH/VL"):
            self.build(contact_cutoff=7.99)

    def test_numbering_class_count_and_map(self):
        for fmt, update in (("VH", {"domain_count": 2}), ("VL", {"chain": "H"}),
                            ("VL", {"domain_count": 0}), ("VH", {"start": -1}),
                            ("VH", {"sequence_indices": {"41": 1, "42": 1, "43": 3}})):
            def bad(sequence, domain_format):
                result = numbered(sequence, domain_format)
                if domain_format == fmt:
                    result.update(update)
                return result
            self.number.side_effect = bad
            with self.subTest(fmt=fmt, update=update), self.assertRaises(ValueError):
                self.build()

    def test_provenance_values_and_copy(self):
        for key, value in (("source_id", ""), ("source_id", 5), ("structure_kind", "validated"),
                           ("assembly_context", "inferred"), ("reviewed", False), ("reviewed", 1)):
            provenance = dict(self.provenance, **{key: value})
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                self.build(provenance=provenance)
        for key in self.provenance:
            provenance = dict(self.provenance)
            del provenance[key]
            with self.assertRaises(ValueError):
                self.build(provenance=provenance)
        self.provenance["notes"] = ["caller assertion"]
        result = self.build()
        self.provenance["notes"].append("changed")
        self.assertEqual(result["provenance"]["notes"], ["caller assertion"])
        with self.assertRaises(ValueError):
            self.build(provenance=dict(self.provenance, score=float("nan")))

    def test_input_limits_cutoffs_and_format(self):
        for sequence in ("", "acdt", "AX", "A" * 2001):
            with self.assertRaises(ValueError):
                module.build_interface_evidence(sequence, "x", vh_chain_id="h", vl_chain_id="l", provenance=self.provenance)
        self.assertEqual(self.number.call_count, 0)
        for cutoff in (0, -1, True, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                self.build(contact_cutoff=cutoff)
            with self.assertRaises(ValueError):
                self.build(antigen_cutoff=cutoff)
        self.build(structure_format="mmcif")
        self.parse.assert_called_with("synthetic structure\n", format="mmcif")
        with self.assertRaises(ValueError):
            self.build(structure_format="auto")
        self.chains[1] = chain("12E8_L", "A" * 2001, [(0, 0, 0)] * 2001)
        with self.assertRaisesRegex(ValueError, "1..2000"):
            self.build()

    def test_native_exceptions_propagate(self):
        self.parse.side_effect = RuntimeError("native failure")
        with self.assertRaisesRegex(RuntimeError, "native failure"):
            self.build()


if __name__ == "__main__":
    unittest.main()
