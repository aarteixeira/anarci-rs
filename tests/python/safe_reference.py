"""Release inputs: public 12e8 anti-lysozyme antibody and synthetic derivatives.

Only the two 12e8 entries from ANARCI 1.3's anarci.py example are retained.
No repository biological datasets are imported by these tests.
"""
from pathlib import Path
import hashlib
import json

SCHEMES = ("imgt", "chothia", "kabat", "martin", "aho", "wolfguy")
SEQUENCES = [
    ("12e8:H", "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAAKTTPPSVYPLAP"),
    ("12e8:L", "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASV"),
]
HEAVY, LIGHT = (sequence for _, sequence in SEQUENCES)
SYNTHETIC = [
    ("synthetic-no-hit", "ACDEFGHIKLMNPQRSTVWY" * 5),
    ("synthetic-truncated-12e8", HEAVY[15:105]),
    ("synthetic-tandem-12e8", HEAVY + "GGGGSGGGGS" + LIGHT),
]


def dataset(count=16, unique=True):
    """Full-length synthetic variants; mutations stay outside the variable domain."""
    if count < 1:
        raise ValueError("count must be positive")
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    result = []
    for index in range(count):
        sequence = SEQUENCES[index % 2][1]
        if unique:
            value, suffix = index, ""
            for _ in range(5):
                suffix += alphabet[value % len(alphabet)]
                value //= len(alphabet)
            sequence += suffix
        result.append(("synthetic-%d" % index, sequence))
    return result


def dataset_sha256(records):
    return hashlib.sha256(json.dumps(records, separators=(",", ":")).encode()).hexdigest()


def expected():
    return json.loads(Path(__file__).with_name("safe_expected.json").read_text())
