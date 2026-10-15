"""Original ANARCI function signatures, implemented with native batch computation."""
from ._api import ALLOW, SPECIES, SCHEMES, _batch, _sequences, _scheme, _hmmerpath, _emit
from . import _native
from .io import *
from .schemes import *
from .schemes import _states

scheme_short_to_long = SCHEMES.copy()
scheme_names = list(scheme_short_to_long)
all_reference_states = list(range(1, 129))
chain_type_to_class = {"H":"H", "K":"L", "L":"L", "A":"A", "B":"B", "G":"G", "D":"D"}


class HMMscanError(Exception):
    """Failure while scanning a receptor sequence with HMMER."""


def anarci(sequences, scheme="imgt", database="ALL", output=False, outfile=None,
           csv=False, allow={"H", "K", "L", "A", "B", "G", "D"}, hmmerpath="",
           ncpu=None, assign_germline=False, allowed_species=["human", "mouse"], bit_score_threshold=80):
    """Number sequence pairs with the reference ANARCI defaults and result shapes."""
    if assign_germline and allowed_species is not None and len(allowed_species) == 0:
        sequences = list(sequences)
        result = _batch(sequences, scheme, database, False, outfile, csv, allow, hmmerpath,
                        ncpu, False, allowed_species, bit_score_threshold)
        if any(result[0]):
            raise ValueError("max() arg is an empty sequence")
        _emit(sequences, result, output, outfile, csv)
        return result
    return _batch(sequences, scheme, database, output, outfile, csv, allow, hmmerpath,
                  ncpu, assign_germline, allowed_species, bit_score_threshold)


def run_anarci(seq, ncpu=1, **kwargs):
    """Read a FASTA, protein string or sequence pairs and number them in Rust."""
    sequences = _sequences(seq, legacy=True)
    return (sequences,) + anarci(sequences, ncpu=ncpu, **kwargs)


def number(sequence, scheme="imgt", database="ALL", allow={"H", "K", "L", "A", "B", "G", "D"},
           allowed_species=["human", "mouse"]):
    """Return first-domain numbering and chain class, or (False, False)."""
    validate_sequence(sequence)
    scheme = _scheme(scheme, number=True)
    if len(sequence) < 70:
        return False, False
    try:
        # The legacy wrapper uses its human/mouse defaults for any false value.
        species = allowed_species if allowed_species else SPECIES
        numbered, details, _ = anarci([("sequence_0", sequence)], scheme=scheme,
                                      database=database, allow=allow, allowed_species=species)
    except AssertionError:
        return False, False
    if numbered[0]:
        return numbered[0][0][0], chain_type_to_class[details[0][0]["chain_type"]]
    return False, False


def run_hmmer(sequence_list, hmm_database="ALL", hmmerpath="", ncpu=None,
              bit_score_threshold=80, hmmer_species=None):
    """Return hit tables, IMGT alignment state vectors and domain details."""
    _hmmerpath(hmmerpath)
    return _native._run_hmmer([tuple(record) for record in sequence_list], str(hmm_database), ncpu, bit_score_threshold, hmmer_species)


def number_sequence_from_alignment(state_vector, sequence, scheme="imgt", chain_type=None):
    """Convert one alignment state vector into a requested numbering scheme."""
    return _native._number_alignment(_states(state_vector), "".join(sequence), scheme, chain_type)


def number_sequences_from_alignment(sequences, alignments, scheme="imgt",
                                    allow={"H", "K", "L", "A", "B", "G", "D"},
                                    assign_germline=False, allowed_species=None):
    """Number supplied alignments without repeating the HMM search."""
    if len(sequences) != len(alignments):
        raise ValueError("Each sequence requires one alignment result")
    numbered, alignment_details, hit_tables = [], [], []
    for name_sequence, (hits, vectors, details) in zip(sequences, alignments):
        if len(vectors) != len(details):
            raise ValueError("Each state vector requires one domain detail record")
        name, sequence = name_sequence
        nums, annotations = [], []
        for vector, detail in zip(vectors, details):
            detail["scheme"], detail["query_name"] = scheme, name
            if vector and detail["chain_type"] in allow:
                nums.append(validate_numbering(number_sequence_from_alignment(vector, sequence, scheme,
                                                                                detail["chain_type"]), name_sequence))
                if assign_germline:
                    detail["germlines"] = run_germline_assignment(vector, sequence, detail["chain_type"], allowed_species)
                annotations.append(detail)
        numbered.append(nums or None)
        alignment_details.append(annotations or None)
        hit_tables.append(hits)
    return numbered, alignment_details, hit_tables


def get_identity(state_sequence, germline_sequence):
    return _native._get_identity("".join(state_sequence), "".join(germline_sequence))


def get_hmm_length(species, ctype):
    return _native._get_hmm_length(species, ctype)


def run_germline_assignment(state_vector, sequence, chain_type, allowed_species=None):
    return _native._germline_assignment(_states(state_vector), "".join(sequence), chain_type, allowed_species)


def check_for_j(sequences, alignments, scheme):
    """Rescan J regions for long CDR3 alignments and update alignments in place."""
    if len(sequences) != len(alignments):
        raise ValueError("Each sequence requires one alignment result")
    for (name, sequence), alignment in zip(sequences, alignments):
        _, vectors, details = alignment
        updated, ends = _native._check_for_j(name, sequence, [_states(vector) for vector in vectors], [d["query_end"] for d in details])
        vectors[:] = updated
        for detail, end in zip(details, ends):
            detail["query_end"] = end


def _domains_are_same(dom1, dom2):
    dom1, dom2 = sorted([dom1, dom2], key=lambda domain: domain.query_start)
    return dom2.query_start < dom1.query_end


def _hsp_fields(hsp):
    return (hsp.hit_id, hsp.hit_description, hsp.evalue, hsp.bitscore, hsp.bias,
            hsp.query_start, hsp.query_end, hsp.hit_start, hsp.hit_end,
            hsp.aln_annotation["RF"], hsp.aln_annotation["PP"], getattr(hsp, "order", 0))


def _parse_hmmer_query(query, bit_score_threshold=80, hmmer_species=None):
    return _native._parse_hsps([_hsp_fields(hsp) for hsp in query.hsps], query.seq_len,
                                bit_score_threshold, hmmer_species)


def _hmm_alignment_to_states(hsp, n, seq_length):
    return _native._hmm_states(_hsp_fields(hsp), n, seq_length)


def parse_hmmer_output(filedescriptor="", bit_score_threshold=80, hmmer_species=None):
    """Parse saved HMMER text with the optional Biopython parser.

    In-process ``run_hmmer`` needs no Biopython installation.
    """
    import os
    try:
        from Bio.SearchIO.HmmerIO import Hmmer3TextParser
    except ImportError:
        raise ImportError("Parsing saved HMMER text requires the optional biopython package") from None
    opener = os.fdopen if isinstance(filedescriptor, int) else open
    with opener(filedescriptor) as handle:
        return [_parse_hmmer_query(query, bit_score_threshold, hmmer_species)
                for query in Hmmer3TextParser(handle)]
