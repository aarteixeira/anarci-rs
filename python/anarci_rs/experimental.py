"""Research codec, augmentation, and dataset helpers backed by Rust.

Labels are ``("match", position)`` or ``("insert", position)`` for positions
1..128, or ``"out"``, ``"linker"``, ``"constant"``, ``"hinge"``, ``"other"``.
States are ``(position, "M" | "I" | "D", sequence_index | None)`` tuples.
Sequence indices are zero based. Sequences use single-byte ASCII characters.

These utilities operate on supplied labels and scores. They do not load a
trained numbering model. Model construction and training require a source
build with the optional ``learned-model`` Cargo feature and run on the CPU.
"""
from dataclasses import dataclass
from . import _native

NUM_LABELS = 261
FULL_THR = 0.95
CDR3_THR = 0.80
MODEL_AVAILABLE = hasattr(_native, "ExperimentalModel")
_OUTSIDE = ("out", "linker", "constant", "hinge", "other")


def label_to_index(label):
    """Return the model vocabulary index for a label; reject invalid positions."""
    if isinstance(label, str) and label in _OUTSIDE:
        value = 256 + _OUTSIDE.index(label)
    elif (isinstance(label, (tuple, list)) and len(label) == 2
          and label[0] in ("match", "insert")
          and isinstance(label[1], int) and not isinstance(label[1], bool)
          and 1 <= label[1] <= 128):
        value = label[1] - 1 + (128 if label[0] == "insert" else 0)
    else:
        raise ValueError("expected ('match'|'insert', 1..128) or a non-domain label")
    return _native.exp_label_to_index(value)


def index_to_label(index):
    """Return a label for a vocabulary index in ``range(NUM_LABELS)``."""
    if not isinstance(index, int) or isinstance(index, bool) or not 0 <= index < NUM_LABELS:
        raise ValueError("label index must be an integer in 0..261")
    index = _native.exp_index_to_label(index)
    if index < 128:
        return ("match", index + 1)
    if index < 256:
        return ("insert", index - 127)
    return _OUTSIDE[index - 256]


def encode_domain(states):
    """Convert states to ``(sequence_index, label)`` pairs, dropping deletes."""
    return [(si, index_to_label(label)) for si, label in _native.exp_encode_domain(list(states))]


def decode_domain(residues):
    """Reconstruct one domain, including deletes between occupied positions.

    Residues are ``(sequence_index, label)`` pairs with increasing indices and
    match positions. Non-domain labels are ignored, as in the Rust codec.
    """
    return _native.exp_decode_domain([(si, label_to_index(label)) for si, label in residues])


def strip_boundary_deletes(states):
    """Remove leading and trailing delete states."""
    return _native.exp_strip_boundary_deletes(list(states))


def decode_sequence(labels):
    """Segment labels into domain state vectors using separators and match resets.

    Abutting partials whose positions keep increasing can merge. Use
    :func:`concat_domains` when a guaranteed separator is needed.
    """
    return _native.exp_decode_sequence([label_to_index(label) for label in labels])


def viterbi_decode(logits):
    """Decode finite ``length x 261`` scores under the Rust transition rules."""
    return [index_to_label(i) for i in _native.exp_viterbi_decode([list(row) for row in logits])]


@dataclass(frozen=True)
class LabeledSeq:
    """ASCII sequence with one label per residue."""
    sequence: str
    labels: tuple

    def __post_init__(self):
        if not isinstance(self.sequence, str) or not self.sequence.isascii():
            raise ValueError("sequence must be an ASCII string")
        labels = tuple(index_to_label(label_to_index(label)) for label in self.labels)
        if len(self.sequence) != len(labels):
            raise ValueError("sequence and labels must have equal lengths")
        object.__setattr__(self, "labels", labels)

    @classmethod
    def from_domain(cls, sequence, states):
        """Extract residues indexed by a single domain state vector."""
        return _segment(_native.exp_from_domain(sequence, list(states)))

    @classmethod
    def region(cls, sequence, label):
        """Label a segment with one non-domain label."""
        return _segment(_native.exp_region(sequence, label_to_index(label)))

    @classmethod
    def junk(cls, sequence):
        """Label every residue as ``out``."""
        return _segment(_native.exp_junk(sequence))

    def restrict_positions(self, lo, hi):
        """Retain domain residues with positions in the inclusive range."""
        return _segment(_native.exp_restrict_positions(_segment_row(self), lo, hi))

    def restrict_window(self, start, end):
        """Retain ``[start, end)``; bounds are clamped to the sequence length."""
        return _segment(_native.exp_restrict_window(_segment_row(self), start, end))

    def domains(self):
        """Return reconstructed state vectors in sequence order."""
        return _native.exp_domains(_segment_row(self))


def _segment(row):
    sequence, labels = row
    return LabeledSeq(sequence, tuple(index_to_label(i) for i in labels))


def _segment_row(value):
    if not isinstance(value, LabeledSeq):
        raise TypeError("expected LabeledSeq")
    return (value.sequence, [label_to_index(label) for label in value.labels])


def concat(segments):
    """Concatenate segments directly, preserving every supplied label."""
    return _segment(_native.exp_concat([_segment_row(s) for s in segments]))


def concat_domains(segments):
    """Concatenate with a ``G`` residue labeled ``out`` between segments."""
    return _segment(_native.exp_concat_domains([_segment_row(s) for s in segments]))


def gly_ser_linker(repeats):
    """Return ``GGGGS`` repeats labeled as linker."""
    return _segment(_native.exp_gly_ser_linker(repeats))


@dataclass(frozen=True)
class DomainRecord:
    """A full source sequence and one domain's state vector.

    Domains with the same full sequence are treated as one source record by
    the clustering and split helpers.
    """
    sequence: str
    states: tuple

    def __post_init__(self):
        # Reuse native validation, including state indices and single-domain IDs.
        states = tuple(tuple(s) for s in self.states)
        _native.exp_from_domain(self.sequence, list(states))
        object.__setattr__(self, "states", states)


def _record_row(record):
    if not isinstance(record, DomainRecord):
        raise TypeError("expected DomainRecord")
    return (record.sequence, list(record.states))


def cdr3_of(record):
    """Extract residues at IMGT positions 105..117 from supplied states."""
    return _native.exp_cdr3_of(_record_row(record))


def overlap_identity(a, b):
    """Return best gap-free overlap identity, normalized by the shorter length."""
    return _native.exp_overlap_identity(a, b)


def cluster_records(records, full_thr=FULL_THR, cdr3_thr=CDR3_THR):
    """Return cluster IDs from full-sequence OR CDR3 overlap similarity."""
    return _native.exp_cluster_records([_record_row(r) for r in records], full_thr, cdr3_thr)


def split_for_rep(representative):
    """Return deterministic ``train``, ``val``, or ``test`` from ASCII bytes."""
    return _native.exp_split_for_rep(representative)


def splits_of(records):
    """Assign clustered source records to deterministic splits at default thresholds."""
    return _native.exp_splits_of([_record_row(r) for r in records])


def near_duplicate(a, b):
    """Test full-sequence or CDR3 similarity at the dataset's default thresholds."""
    return _native.exp_near_duplicate(_record_row(a), _record_row(b))


def build_dataset(records, partial_windows=()):
    """Build labeled single, partial, multidomain, and built-in junk examples.

    Each returned dict contains ``data`` (LabeledSeq), ``split``, ``clusters``,
    and ``kind``. This calls the existing Rust research builder, including its
    fixed junk negatives and deterministic clustering/split rules.
    """
    rows = _native.exp_build_dataset([_record_row(r) for r in records], list(partial_windows))
    return [dict(data=_segment(data), split=split, clusters=clusters, kind=kind)
            for data, split, clusters, kind in rows]


def _require_model():
    if not MODEL_AVAILABLE:
        raise RuntimeError(
            "CPU research model tools require a source build with "
            "`maturin build --release --features learned-model`; "
            "the default wheel includes the codec and dataset tools."
        )


class Model:
    """Untrained CPU transformer; optional ``learned-model`` source build.

    ``dim`` and ``ff`` must be positive; ``layers`` can be zero. Input sequences
    have 1..1024 ASCII residues. No trained weights are loaded by this class.
    """

    def __init__(self, dim=64, layers=2, ff=128):
        _require_model()
        self._inner = _native.ExperimentalModel(dim, layers, ff)

    @classmethod
    def new(cls, dim=64, layers=2, ff=128):
        """Construct an untrained model with the supplied dimensions."""
        return cls(dim=dim, layers=layers, ff=ff)

    @classmethod
    def _from_native(cls, inner):
        instance = cls.__new__(cls)
        instance._inner = inner
        return instance

    def logits(self, sequence):
        """Return per-residue scores shaped ``[length, 261]``."""
        return self._inner.logits(sequence)

    def forward(self, tokens):
        """Return ``[batch, length, 261]`` scores for a rectangular token batch.

        Token IDs are 0..21, with 0 used for padding. The existing research
        model attends to all supplied positions, including padding tokens.
        """
        return self._inner.forward([list(row) for row in tokens])


def tokenize(sequence):
    """Map ASCII residues to token IDs 1..20 and unknown residues to 21."""
    _require_model()
    return _native.exp_tokenize(sequence)


def overfit(examples, dim=64, layers=2, ff=128, steps=1, lr=0.001):
    """Train on supplied examples and return ``(Model, training_accuracy)``.

    Examples are ``(sequence, labels)`` pairs or LabeledSeq instances. Labels
    use the same tuple/string representation as the codec. The CPU optimizer
    is the existing Rust AdamW helper. Accuracy describes this training set.
    ``steps=0`` returns the initial model's training accuracy.
    """
    _require_model()
    rows = []
    for example in examples:
        sequence, labels = ((example.sequence, example.labels)
                            if isinstance(example, LabeledSeq) else example)
        rows.append((sequence, [label_to_index(label) for label in labels]))
    inner, accuracy = _native.exp_overfit(rows, dim, layers, ff, steps, lr)
    return Model._from_native(inner), accuracy


__all__ = [
    "NUM_LABELS", "FULL_THR", "CDR3_THR", "LabeledSeq", "DomainRecord",
    "label_to_index", "index_to_label", "encode_domain", "decode_domain",
    "strip_boundary_deletes", "decode_sequence", "viterbi_decode", "concat",
    "concat_domains", "gly_ser_linker", "cdr3_of", "overlap_identity",
    "cluster_records", "split_for_rep", "splits_of", "near_duplicate", "build_dataset",
    "MODEL_AVAILABLE", "Model", "tokenize", "overfit",
]
