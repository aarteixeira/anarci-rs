"""Bundled ANARCI germline reference data, with original dictionary structure."""
import json
from anarci_rs import _native
_data = json.loads(_native._germlines_json())
all_germlines = {
    segment: {chain: {species: dict(genes) for species, genes in species_map.items()}
              for chain, species_map in chains.items()}
    for segment, chains in _data["germlines"].items()
}
del _data
all_species = list(all_germlines["V"]["H"])
