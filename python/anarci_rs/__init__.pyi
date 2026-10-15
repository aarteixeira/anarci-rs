from typing import Any, Iterable, Optional, Sequence, TextIO, Tuple, List, Dict, Union
from os import PathLike

__version__: str
all_germlines: Dict[str, Any]
all_species: List[str]

def anarci(sequences: Iterable[Tuple[str, str]], scheme: str = ..., database: Union[str, PathLike] = ...,
           output: bool = ..., outfile: Optional[Union[str, PathLike]] = ..., csv: bool = ...,
           allow: Optional[Iterable[str]] = ..., hmmerpath: str = ..., ncpu: Optional[int] = ...,
           assign_germline: bool = ..., allowed_species: Optional[Sequence[str]] = ...,
           bit_score_threshold: float = ..., germline_method: Optional[str] = ...,
           annotate_regions: bool = ...) -> Tuple[List[Any], List[Any], List[Any]]: ...
def run_anarci(seq: Any, scheme: str = ..., database: Union[str, PathLike] = ...,
               output: bool = ..., outfile: Optional[Union[str, PathLike]] = ..., csv: bool = ...,
               allow: Optional[Iterable[str]] = ..., hmmerpath: str = ..., ncpu: int = ...,
               assign_germline: bool = ..., allowed_species: Optional[Sequence[str]] = ...,
               bit_score_threshold: float = ..., germline_method: Optional[str] = ...,
               annotate_regions: bool = ...) -> Tuple[Any, List[Any], List[Any], List[Any]]: ...
def number(sequence: str, scheme: str = ..., database: Union[str, PathLike] = ...,
           allow: Optional[Iterable[str]] = ..., allowed_species: Optional[Sequence[str]] = ...,
           min_length: int = ..., bit_score_threshold: float = ...) -> Tuple[Any, Any]: ...
from .io import read_fasta as read_fasta, fasta_iter as fasta_iter, write_fasta as write_fasta
from .io import validate_sequence as validate_sequence, validate_numbering as validate_numbering
from .io import anarci_output as anarci_output, csv_output as csv_output
from .compat import run_hmmer as run_hmmer, number_sequence_from_alignment as number_sequence_from_alignment
from .compat import number_sequences_from_alignment as number_sequences_from_alignment
from .compat import get_identity as get_identity, get_hmm_length as get_hmm_length
from .compat import run_germline_assignment as run_germline_assignment, check_for_j as check_for_j
from .compat import parse_hmmer_output as parse_hmmer_output
from ._native import read_sequences as read_sequences, translate_dna as translate_dna
from ._native import parse_structure as parse_structure, annotate_protein as annotate_protein
from ._native import annotate_dna as annotate_dna, pair_structure as pair_structure
from ._native import format_annotations as format_annotations
