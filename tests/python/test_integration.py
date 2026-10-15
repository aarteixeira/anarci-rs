"""Independent release checks using only the permitted antibody fixture."""
import concurrent.futures
import csv
import importlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from safe_reference import HEAVY, LIGHT, SCHEMES, SEQUENCES, SYNTHETIC, expected


def plain(value):
    return json.loads(json.dumps(value))


def pdb_fixture(antigen=False):
    """Synthetic C-alpha geometry, not an experimental antibody structure."""
    names = dict(zip("ACDEFGHIKLMNPQRSTVWY", "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()))
    lines = []
    chains = [("H", HEAVY, 0), ("L", LIGHT, 5)]
    if antigen:
        chains.append(("A", "ACDEFGHIKLMNPQRSTVWY" * 5, 2.5))
    for chain, sequence, y in chains:
        for i, aa in enumerate(sequence, 1):
            lines.append("ATOM  %5d  CA  %3s %s%4d    %8.3f%8.3f%8.3f  1.00 20.00           C" %
                         (len(lines) + 1, names[aa], chain, i, i * 3.8, y, 0))
    return "\n".join(lines) + "\nEND\n"


def dna_fixture():
    codons = dict(zip("ACDEFGHIKLMNPQRSTVWY", "GCT TGT GAT GAA TTT GGT CAT ATT AAA CTG ATG AAT CCT CAA CGT TCT ACT GTT TGG TAT".split()))
    return "".join(codons[aa] for aa in HEAVY)


def cif_fixture():
    columns = ("group_PDB", "label_atom_id", "label_comp_id", "auth_asym_id", "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z")
    text = "data_synthetic\nloop_\n" + "\n".join("_atom_site." + column for column in columns) + "\n"
    for line in pdb_fixture().splitlines():
        if line.startswith("ATOM"):
            text += "ATOM CA %s %s %s %s %s %s\n" % tuple(line[a:b].strip() for a, b in (
                (17, 20), (21, 22), (22, 26), (30, 38), (38, 46), (46, 54)))
    return text + "#\n"


class ReferenceIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.extended = importlib.import_module("anarci_rs")
        cls.compat = importlib.import_module("anarci")
        cls.oracle = expected()

    def test_scheme_helpers_against_pinned_numbering(self):
        schemes = importlib.import_module("anarci.schemes")
        api = importlib.import_module("anarci.anarci")
        alignments = json.loads(Path(__file__).with_name("safe_alignment.json").read_text())
        for index, (_, sequence) in enumerate(SEQUENCES):
            vector = [(tuple(state), query_index) for state, query_index in alignments[index][1][0]]
            alignments[index][1][0] = vector
            germline = api.run_germline_assignment(vector, sequence, "H" if index == 0 else "K", ["human", "mouse"])
            self.assertEqual(plain(germline), self.oracle["imgt"][1][index][0]["germlines"])
            for scheme in SCHEMES:
                reference = self.oracle[scheme][0][index][0]
                direct = api.number_sequence_from_alignment(vector, sequence, scheme=scheme, chain_type="H" if index == 0 else "K")
                self.assertEqual(plain(direct), reference)
                if scheme == "imgt":
                    result = schemes.number_imgt(vector, sequence)
                elif scheme == "aho":
                    result = schemes.number_aho(vector, sequence, "H" if index == 0 else "K")
                else:
                    function = getattr(schemes, "number_%s_%s" % (scheme, "heavy" if index == 0 else "light"))
                    result = function(vector, sequence)
                self.assertEqual(plain(result), reference)
        api.check_for_j(SEQUENCES, alignments, "imgt")
        numbered = api.number_sequences_from_alignment(SEQUENCES, alignments, scheme="imgt",
                                                        assign_germline=True, allowed_species=["human", "mouse"])
        self.assertEqual(plain(numbered[0]), self.oracle["imgt"][0])

    def test_numeric_helpers_against_pinned_reference(self):
        cases = json.loads(Path(__file__).with_name("safe_helpers.json").read_text())
        for case in cases:
            module = importlib.import_module("anarci." + case["module"])
            arguments = case["args"]
            if case["name"] == "smooth_insertions":
                arguments[0] = [(tuple(state), query_index) for state, query_index in arguments[0]]
            if case["name"] == "_get_wolfguy_L1":
                arguments[0] = [(tuple(position), aa) for position, aa in arguments[0]]
            with self.subTest(helper=case["name"], arguments=arguments[1:]):
                self.assertEqual(plain(getattr(module, case["name"])(*arguments)), case["expected"])

    @unittest.skipUnless(importlib.util.find_spec("Bio"), "optional biopython parser requires biopython")
    def test_saved_hmmer_text_parser_against_reference_alignment(self):
        api = importlib.import_module("anarci.anarci")
        fixture = Path(__file__).with_name("safe_hmmer.txt")
        reference = json.loads(fixture.with_name("safe_alignment.json").read_text())[0]
        parsed = api.parse_hmmer_output(str(fixture))
        self.assertEqual(plain(parsed), [reference])

    def test_exact_all_schemes_numbering_identity_and_evalues(self):
        for scheme in SCHEMES:
            with self.subTest(scheme=scheme):
                actual = plain(self.compat.anarci(SEQUENCES, scheme=scheme, assign_germline=True, ncpu=1))
                expected_result = self.oracle[scheme]
                self.assertEqual(actual[0], expected_result[0])
                for actual_domains, expected_domains in zip(actual[1], expected_result[1]):
                    self.assertEqual(len(actual_domains), len(expected_domains))
                    for detail, reference in zip(actual_domains, expected_domains):
                        for key in ("id", "description", "species", "chain_type", "scheme", "query_name", "query_start", "query_end"):
                            self.assertEqual(detail[key], reference[key], key)
                        self.assertTrue(math.isclose(float(detail["evalue"]), float(reference["evalue"]), rel_tol=0.02, abs_tol=0.0))
                        for key in ("bitscore", "bias"):
                            self.assertAlmostEqual(detail[key], reference[key], delta=0.051)
                        for gene in ("v_gene", "j_gene"):
                            self.assertEqual(detail["germlines"][gene][0], reference["germlines"][gene][0])
                            self.assertAlmostEqual(detail["germlines"][gene][1], reference["germlines"][gene][1], delta=1e-12)
                self.assertEqual(actual[2][0][0], expected_result[2][0][0])
                self.assertEqual(actual[2][0][1][0], expected_result[2][0][1][0])

    def test_engine_switching_methods_filters_and_synthetic_domains(self):
        for database in ("pan", "ALL", "pan", "ALL"):
            for scheme in SCHEMES:
                for method in ("identity", "evalue"):
                    with self.subTest(database=database, scheme=scheme, method=method):
                        result = self.extended.anarci(SEQUENCES, scheme=scheme, database=database,
                                                     assign_germline=True, germline_method=method,
                                                     allowed_species=["mouse"], ncpu=2)
                        self.assertTrue(all(result[0]))
                        self.assertEqual([d[0]["chain_type"] for d in result[1]], ["H", "K"])
                        self.assertTrue(all(d[0]["germlines"] for d in result[1]))
            filtered = self.extended.anarci(SEQUENCES, database=database, allow={"H"})
            self.assertTrue(filtered[0][0])
            self.assertIsNone(filtered[0][1])
        result = self.extended.anarci(SYNTHETIC, database="ALL", bit_score_threshold=20)
        self.assertIsNone(result[0][0])
        self.assertTrue(result[0][1])
        self.assertEqual(len(result[0][2]), 2)

    def test_seven_chain_classes_with_model_derived_synthetic_inputs(self):
        fixture = json.loads(Path(__file__).with_name("safe_consensus.json").read_text())
        records = [tuple(record) for record in fixture["records"]]
        oracle = fixture["expected_imgt"]
        for database in ("ALL", "pan"):
            result = self.extended.anarci(records, database=database, allowed_species=["human"],
                                         assign_germline=True, germline_method="identity", ncpu=2)
            self.assertEqual([details[0]["chain_type"] for details in result[1]], list("HKLABGD"))
            if database == "ALL":
                self.assertEqual(plain(result[0]), oracle[0])
            for index, chain in enumerate("HKLABGD"):
                filtered = self.extended.anarci([records[index]], database=database, allow={chain})
                self.assertTrue(filtered[0][0])
        # AHo and IMGT support all seven classes; antibody schemes accept H/K/L.
        for scheme in ("aho", "chothia", "kabat", "martin", "wolfguy"):
            selected = records if scheme == "aho" else records[:3]
            result = self.extended.anarci(selected, database="ALL", scheme=scheme)
            self.assertTrue(all(result[0]))

    def test_concurrent_first_use_in_clean_process(self):
        script = '''
import concurrent.futures
import anarci_rs
sequence = %r
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    values = list(pool.map(lambda i: anarci_rs.number(sequence, database="ALL" if i %% 2 else "pan"), range(12)))
assert all(value[1] == "H" for value in values)
''' % HEAVY
        run = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, timeout=60)
        self.assertEqual(run.returncode, 0, run.stderr)

    def test_negative_inputs_cannot_abort_interpreter(self):
        # A separate interpreter makes a native abort an ordinary regression failure.
        cases = ["'A\\0A'", "'A' * 1000001", "'ACD!'", "'ACD☃'"]
        for expression in cases:
            with self.subTest(expression=expression):
                script = '''
import anarci_rs
try:
    anarci_rs.anarci([("invalid", %s)])
except (ValueError, AssertionError, TypeError):
    pass
else:
    raise RuntimeError("invalid input accepted")
print("rejected")
''' % expression
                run = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, timeout=30)
                self.assertEqual(run.returncode, 0, run.stderr)
                self.assertEqual(run.stdout.strip(), "rejected")

    def test_native_sequence_translation_annotation_and_output(self):
        api = self.extended
        fasta = ">safe description\n" + HEAVY + "\n"
        self.assertEqual(api.read_sequences(fasta)[0]["sequence"], HEAVY)
        fastq = "@safe description\n" + HEAVY + "\n+\n" + "I" * len(HEAVY) + "\n"
        self.assertEqual(api.read_sequences(fastq)[0]["sequence"], HEAVY)
        for malformed in ("@x\nAAAA\n-\nIIII\n", "@x\nAAAA\n+\nIII\n"):
            with self.assertRaises((ValueError, AssertionError)):
                api.read_sequences(malformed)
        dna = dna_fixture()
        self.assertEqual(api.translate_dna(dna, frame=0), HEAVY)
        for frame in range(6):
            self.assertIsInstance(api.translate_dna(dna, frame=frame), str)
        for database in ("pan", "ALL"):
            protein = api.annotate_protein(SEQUENCES, database=database, assign_germline=True)
            self.assertEqual(len(protein), 2)
            self.assertTrue(protein[0]["numbering"])
            dna_result = api.annotate_dna([("synthetic-dna", dna)], database=database)
            self.assertEqual(dna_result["domains"][0]["frame"], 0)
            self.assertEqual(sum(dna_result["frame_counts"]), 1)
        for output_format in ("json", "tsv", "anarci"):
            output = api.format_annotations(SEQUENCES, output_format=output_format)
            self.assertTrue(output)
            if output_format == "json":
                self.assertEqual(len(json.loads(output)), 2)
        identifier = 'safe\ttwo "quotes"\nsecond line'
        rendered = api.format_annotations([(identifier, HEAVY)], output_format="tsv")
        rows = list(csv.reader(io.StringIO(rendered), delimiter="\t"))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[1][0], identifier)
        self.assertEqual(len(rows[0]), len(rows[1]))
        with self.assertRaises((ValueError, AssertionError)):
            api.translate_dna("ATG", frame=6)
        with self.assertRaises((ValueError, AssertionError)):
            api.annotate_dna([("invalid", "ATG!")])

    def test_synthetic_structure_parsing_pairing_and_contacts(self):
        api = self.extended
        text = pdb_fixture()
        chains = api.parse_structure(text, format="pdb")
        self.assertEqual([chain["sequence"] for chain in chains], [HEAVY, LIGHT])
        contacts = []
        for scheme in ("imgt", "kabat"):
            pairing = api.pair_structure(pdb_fixture(antigen=True), format="pdb", scheme=scheme)
            self.assertEqual(pairing["n_domains"], 2)
            fvs = [fv for assembly in pairing["assemblies"] for fv in assembly["fvs"]]
            self.assertEqual(len(fvs), 1)
            self.assertGreater(fvs[0]["contacts"], 8)
            synthetic_partner = [partner for partner in fvs[0]["antigen"] if partner["chain"].endswith("_A")]
            self.assertEqual(len(synthetic_partner), 1)
            self.assertFalse(synthetic_partner[0]["is_antibody"])
            contacts.append(synthetic_partner[0]["contacts"])
        self.assertEqual(contacts[0], contacts[1])
        cif = "data_synthetic\nloop_\n" + "\n".join("_atom_site." + c for c in (
            "group_PDB", "label_atom_id", "label_comp_id", "auth_asym_id", "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z"))
        cif += "\nATOM CA GLU H 1 1 2 3\nATOM CA VAL H 2 4 5 6\n#\n"
        self.assertEqual(api.parse_structure(cif, format="mmcif")[0]["sequence"], "EV")
        wrapped = cif.replace("CA GLU H 1 1 2 3", "'CA' GLU H\n1 1 2 3")
        self.assertEqual(api.parse_structure(wrapped, format="cif")[0]["sequence"], "EV")
        missing_coords = text.splitlines()[0]
        missing_coords = missing_coords[:30] + " " * 24 + missing_coords[54:] + "\n"
        self.assertEqual(api.parse_structure(missing_coords)[0]["residues"][0]["ca"], [None, None, None])
        with self.assertRaises((ValueError, AssertionError)):
            api.parse_structure("", format="unsupported")


class InstalledCli(unittest.TestCase):
    def cli(self, *arguments, text=None, ok=True):
        executable = shutil.which("anarci-rs", path=str(Path(sys.executable).parent))
        self.assertIsNotNone(executable, "installed wheel must supply anarci-rs")
        run = subprocess.run([executable, *arguments], input=text, capture_output=True, text=True, timeout=60)
        if ok:
            self.assertEqual(run.returncode, 0, run.stderr)
        else:
            self.assertGreater(run.returncode, 0, run.stderr)
        return run

    def test_help_stdin_output_errors_and_precedence(self):
        self.assertIn("--input-type", self.cli("--help").stdout)
        fasta = ">safe\n" + HEAVY + "\n"
        rows = json.loads(self.cli("--input-type", "protein", "--format", "json", text=fasta).stdout)
        self.assertEqual(rows[0]["chain_type"], "H")
        self.cli("--unknown-release-test-flag", ok=False)
        self.cli("--scheme", "invalid", text=fasta, ok=False)
        self.cli("--input-type", "dna", text=fasta, ok=False)
        self.cli("--pair", "--input-type", "protein", text=fasta, ok=False)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            # An explicit input type takes precedence over a structure filename suffix.
            source = path / "synthetic.pdb"
            source.write_text(fasta)
            output = path / "result.json"
            self.cli(str(source), "--input-type", "protein", "--format", "json", "--output", str(output))
            self.assertEqual(len(json.loads(output.read_text())), 1)
            source.write_text(pdb_fixture())
            paired = json.loads(self.cli(str(source), "--pair", "--format", "json").stdout)
            self.assertEqual(paired["n_domains"], 2)
            source.write_text(cif_fixture())
            for input_type in ("cif", "mmcif"):
                paired = json.loads(self.cli(str(source), "--input-type", input_type, "--pair", "--format", "json").stdout)
                self.assertEqual(paired["n_domains"], 2)


if __name__ == "__main__":
    unittest.main()
