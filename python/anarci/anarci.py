"""ANARCI compatibility entry points and helpers."""
import os
from anarci_rs.compat import *
from .germlines import all_germlines

all_species = list(all_germlines["V"]["H"])
anarci_path = os.path.dirname(__file__)
# Databases live in the native extension; use database="ALL" or a custom HMM path.
HMM_path = None
from anarci_rs.compat import _domains_are_same, _parse_hmmer_query, _hmm_alignment_to_states
