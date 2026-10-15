"""ANARCI scheme helpers backed by the Rust numbering implementation."""
from . import _native

alphabet = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ") + [letter * 2 for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"] + [" "]
blosum62 = dict(_native._blosum62())


def _states(state_vector):
    return [(tuple(state), query_index) for state, query_index in state_vector]


def smooth_insertions(state_vector):
    return _native._smooth_insertions(_states(state_vector))


def number_imgt(state_vector, sequence):
    return _native._number_alignment(_states(state_vector), "".join(sequence), "imgt", None)


def number_aho(state_vector, sequence, chain_type):
    return _native._number_alignment(_states(state_vector), "".join(sequence), "aho", chain_type)


def number_chothia_heavy(state_vector, sequence):
    return _native._number_alignment(_states(state_vector), "".join(sequence), "chothia", "H")


def number_chothia_light(state_vector, sequence):
    return _native._number_alignment(_states(state_vector), "".join(sequence), "chothia", "L")


def number_kabat_heavy(state_vector, sequence):
    return _native._number_alignment(_states(state_vector), "".join(sequence), "kabat", "H")


def number_kabat_light(state_vector, sequence):
    return _native._number_alignment(_states(state_vector), "".join(sequence), "kabat", "L")


def number_martin_heavy(state_vector, sequence):
    return _native._number_alignment(_states(state_vector), "".join(sequence), "martin", "H")


def number_martin_light(state_vector, sequence):
    return _native._number_alignment(_states(state_vector), "".join(sequence), "martin", "L")


def number_wolfguy_heavy(state_vector, sequence):
    return _native._number_alignment(_states(state_vector), "".join(sequence), "wolfguy", "H")


def number_wolfguy_light(state_vector, sequence):
    return _native._number_alignment(_states(state_vector), "".join(sequence), "wolfguy", "L")


def get_imgt_cdr(length, maxlength, start, end):
    return _native._get_imgt_cdr(length, maxlength, start, end)


def get_cdr3_annotations(length, scheme="imgt", chain_type=""):
    if scheme == "imgt":
        assert length - 13 < 50, "Too many insertions for numbering scheme to handle"
        # Legacy CDR3 inserts alternate from position 111, unlike get_imgt_cdr.
        annotations = [None] * max(length, 13)
        front, back = 0, -1
        for index in range(min(length, 13)):
            if index % 2:
                annotations[back] = (118 + back, " ")
                back -= 1
            else:
                annotations[front] = (105 + front, " ")
                front += 1
        letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        for index in range(max(0, length - 13)):
            if index % 2:
                annotations[back] = (112, letters[::-1][back + 6])
                back -= 1
            else:
                annotations[front] = (111, letters[front - 7])
                front += 1
        return annotations
    return _native._get_cdr3_annotations(length, scheme, chain_type)


def gap_missing(numbering):
    result = []
    last = 0
    for region in numbering:
        for (position, insertion), amino_acid in region:
            result.extend(((index, " "), "-") for index in range(last + 1, position))
            result.append(((position, insertion), amino_acid))
            last = position
    return result


def _get_wolfguy_L1(seq, length):
    return _native._get_wolfguy_l1([(tuple(position), aa) for position, aa in seq], length)


__all__ = ["alphabet", "blosum62", "smooth_insertions", "number_imgt", "number_aho",
           "number_chothia_heavy", "number_chothia_light", "number_kabat_heavy",
           "number_kabat_light", "number_martin_heavy", "number_martin_light",
           "number_wolfguy_heavy", "number_wolfguy_light", "get_imgt_cdr",
           "get_cdr3_annotations", "gap_missing"]


def _number_regions(sequence, state_vector, state_string, region_string, region_index_dict,
                    rels, n_regions, exclude_deletions):
    return _native._number_regions("".join(sequence), _states(state_vector), state_string, region_string,
                                    region_index_dict, rels, n_regions, exclude_deletions)
