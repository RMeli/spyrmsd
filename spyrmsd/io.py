try:
    from spyrmsd.optional.obabel import (
        load,
        loadall,
        adjacency_matrix,
        to_molecule,
        numatoms,
        numbonds,
        bonds,
    )

except ImportError:
    from spyrmsd.optional.rdkit import (
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
