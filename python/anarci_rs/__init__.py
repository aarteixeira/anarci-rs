"""Antigen receptor numbering with in-process Rust computation.

The extended API defaults to the pan database. Use ``import anarci`` for the
original ANARCI signatures and ALL database defaults.
"""
from importlib.metadata import PackageNotFoundError, version
from . import _native
from ._api import ALLOW, _batch, _sequences
from .io import read_fasta, fasta_iter, write_fasta, validate_sequence, validate_numbering, anarci_output, csv_output

try:
    __version__ = version("anarci-rs")
except PackageNotFoundError:
    __version__ = "0+local"


def anarci(sequences, scheme="imgt", database="pan", output=False, outfile=None,
           csv=False, allow=None, hmmerpath="", ncpu=None, assign_germline=False,
           allowed_species=None, bit_score_threshold=80.0, germline_method=None,
           annotate_regions=False):
    """Number a batch of (name, protein) pairs; return numbering, details and hits.

    ``database`` accepts ``pan``, ``ALL`` or a custom HMM file path. Output can
    target stdout, an ANARCI text file, or per-chain CSV files. Native batch
    computation releases the GIL; ``ncpu`` selects Rust worker threads.
    """
    return _batch(sequences, scheme, database, output, outfile, csv, allow, hmmerpath,
                  ncpu, assign_germline, allowed_species, bit_score_threshold,
                  germline_method, annotate_regions)


def run_anarci(seq, scheme="imgt", database="pan", output=False, outfile=None,
               csv=False, allow=None, hmmerpath="", ncpu=1, assign_germline=False,
               allowed_species=None, bit_score_threshold=80.0, germline_method=None,
               annotate_regions=False):
    """Read protein input and return sequences, numbering, details and hits."""
    sequences = _sequences(seq)
    result = anarci(sequences, scheme, database, output, outfile, csv, allow, hmmerpath,
                    ncpu, assign_germline, allowed_species, bit_score_threshold,
                    germline_method, annotate_regions)
    return (sequences,) + result


def number(sequence, scheme="imgt", database="pan", allow=None,
           allowed_species=["human", "mouse"], min_length=70, bit_score_threshold=80.0):
    """Number the first domain; warn and return (False, False) when none is found.

    ``min_length`` and ``bit_score_threshold`` allow partial-domain searches.
    """
    return _native.number(sequence, scheme, str(database), allow, allowed_species,
                          min_length, bit_score_threshold)


__all__ = ["anarci", "run_anarci", "number", "read_fasta", "fasta_iter", "write_fasta",
           "validate_sequence", "validate_numbering", "anarci_output", "csv_output"]

# Shared utilities keep the same behavior in the extended and compatibility APIs.
from .compat import run_hmmer, number_sequence_from_alignment, number_sequences_from_alignment
from .compat import run_germline_assignment, get_identity, get_hmm_length, check_for_j, parse_hmmer_output
from ._native import read_sequences, translate_dna, parse_structure, annotate_protein, annotate_dna, pair_structure, format_annotations
__all__ += ["run_hmmer", "number_sequence_from_alignment", "number_sequences_from_alignment",
            "run_germline_assignment", "get_identity", "get_hmm_length", "check_for_j", "parse_hmmer_output",
            "read_sequences", "translate_dna", "parse_structure", "annotate_protein", "annotate_dna", "pair_structure", "format_annotations"]
from .germlines import all_germlines, all_species
__all__ += ["all_germlines", "all_species"]
