"""Shared Python boundary for the native batch engine."""
import os
import sys
from . import _native
from .io import anarci_output, csv_output, read_fasta, validate_sequence

ALLOW = {"H", "K", "L", "A", "B", "G", "D"}
SPECIES = ["human", "mouse"]
SCHEMES = {"m":"martin", "c":"chothia", "k":"kabat", "imgt":"imgt", "kabat":"kabat",
           "chothia":"chothia", "martin":"martin", "i":"imgt", "a":"aho", "aho":"aho", "wolfguy":"wolfguy", "w":"wolfguy"}


def _scheme(scheme, number=False):
    try:
        return SCHEMES[scheme.lower()]
    except KeyError:
        wording = "Unrecognised to unimplemented" if number else "Unrecognised or unimplemented"
        raise AssertionError(f"{wording} scheme: {scheme}") from None


def _output_options(csv, outfile):
    if csv:
        assert outfile, "If csv output is True then an outfile must be specified"
        path = os.path.dirname(outfile)
        assert not path or os.path.exists(path), f"Output directory {path} does not exist"


def _hmmerpath(hmmerpath):
    if hmmerpath:
        raise NotImplementedError("hmmerpath selects an external HMMER executable; anarci-rs uses the bundled in-process engine")


def _emit(sequences, result, output, outfile, csv):
    if output:
        numbered, details, _ = result
        if csv:
            csv_output(sequences, numbered, details, os.fspath(outfile))
        elif outfile:
            with open(outfile, "w") as handle:
                anarci_output(numbered, sequences, details, handle)
        else:
            anarci_output(numbered, sequences, details, sys.stdout)


def _sequences(seq, legacy=False):
    if isinstance(seq, (list, tuple)):
        assert all(len(item) == 2 for item in seq), "If list or tuple supplied as input format must be [ ('ID1','seq1'), ('ID2', 'seq2'), ... ]"
        return seq if legacy else [(name, sequence) for name, sequence in seq]
    if isinstance(seq, os.PathLike) or isinstance(seq, str) and os.path.isfile(seq):
        return read_fasta(os.fspath(seq))
    if isinstance(seq, str):
        validate_sequence(seq)
        return [["Input sequence", seq]] if legacy else [("Input sequence", seq)]
    if not legacy:
        return _sequences(list(seq))
    raise TypeError("Input must be a sequence string, FASTA path, or list of (name, sequence) pairs")


def _batch(sequences, scheme, database, output, outfile, csv, allow, hmmerpath,
           ncpu, assign_germline, allowed_species, bit_score_threshold,
           germline_method=None, annotate_regions=False):
    scheme = _scheme(scheme)
    if not isinstance(sequences, (list, tuple)):
        sequences = list(sequences)
    _output_options(csv, outfile)
    _hmmerpath(hmmerpath)
    native_sequences = [tuple(record) for record in sequences]
    result = _native.anarci(native_sequences, scheme=scheme, database=os.fspath(database), allow=allow,
                            ncpu=ncpu, assign_germline=assign_germline,
                            allowed_species=allowed_species, bit_score_threshold=bit_score_threshold,
                            germline_method=germline_method, annotate_regions=annotate_regions)
    _emit(sequences, result, output, outfile, csv)
    return result
