try:
    from pyrmsd.optional.obabel import (
        load,
        loadall,
        adjacency_matrix,
        to_molecule,
        numatoms,
        numbonds,
        bonds,
    )

except ImportError:
    from pyrmsd.optional.rdkit import (
        load,
        loadall,
        adjacency_matrix,
        to_molecule,
        numatoms,
        numbonds,
        bonds,
    )

__all__ = [
    "load",
    "loadall",
    "adjacency_matrix",
    "to_molecule",
    "numatoms",
    "numbonds",
    "bonds",
]
