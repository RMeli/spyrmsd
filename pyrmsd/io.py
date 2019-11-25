try:
    from pyrmsd.optional.obabel import (
        load,
        loadall,
        adjacency_matrix,
        to_molecule,
    )

    __all__ = ["load", "loadall", "adjacency_matrix", "to_molecule"]

except ImportError:

    raise NotImplementedError

    from pyrmsd.optional.rdkit import load
