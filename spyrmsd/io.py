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

    try:
        from spyrmsd.optional.rdkit import (
            load,
            loadall,
            adjacency_matrix,
            to_molecule,
            numatoms,
            numbonds,
            bonds,
        )
    except ImportError:
        # Use sPyRMSD as standalone library
        __all__ = []
    else:
        # Avoid flake8 complaint "imported but unused"
        __all__ = [
            "load",
            "loadall",
            "adjacency_matrix",
            "to_molecule",
            "numatoms",
            "numbonds",
            "bonds",
        ]
else:
    # Avoid flake8 complaint "imported but unused"
    __all__ = [
        "load",
        "loadall",
        "adjacency_matrix",
        "to_molecule",
        "numatoms",
        "numbonds",
        "bonds",
    ]
