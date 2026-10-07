"""
Optional xyzgraph-based adjacency builder.
"""

from typing import List, Sequence, Tuple

import numpy as np

from spyrmsd import constants


def _atomic_number_to_symbol(atomic_number: int) -> str:
    try:
        return constants.anum_to_symbol[atomic_number]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported atomic number for xyzgraph adapter: {atomic_number}"
        ) from exc


def _atoms_from_arrays(
    aprops: Sequence[int], coordinates: np.ndarray
) -> List[Tuple[str, Tuple[float, float, float]]]:
    n_atoms = len(aprops)

    assert coordinates.shape == (n_atoms, 3)

    atoms = []
    for atomic_number, position in zip(aprops, coordinates):
        symbol = _atomic_number_to_symbol(int(atomic_number))
        xyz_position = (
            float(position[0]),
            float(position[1]),
            float(position[2]),
        )
        atoms.append((symbol, xyz_position))

    return atoms


def adjacency_matrix_from_atomic_coordinates(
    aprops: np.ndarray, coordinates: np.ndarray
) -> np.ndarray:
    """
    Compute an adjacency matrix from atomic coordinates using xyzgraph.

    Parameters
    ----------
    aprops: numpy.ndarray
        Atomic numbers.
    coordinates: numpy.ndarray
        Atomic coordinates.
    Returns
    -------
    numpy.ndarray
        Adjacency matrix.
    """

    try:
        from xyzgraph import build_graph
    except ImportError as exc:
        raise ImportError(
            "xyzgraph is required for the xyzgraph graph builder. "
            "Install with `pip install spyrmsd[xyzgraph]`."
        ) from exc

    atoms = _atoms_from_arrays(aprops, coordinates)
    graph = build_graph(atoms)

    n_atoms = len(aprops)
    adjacency = np.zeros((n_atoms, n_atoms), dtype=int)

    for i, j in graph.edges():
        ii = int(i)
        jj = int(j)
        if ii < 0 or jj < 0 or ii >= n_atoms or jj >= n_atoms:
            raise ValueError("xyzgraph adapter received out-of-range node indices.")
        adjacency[ii, jj] = adjacency[jj, ii] = 1

    return adjacency


__all__ = ["adjacency_matrix_from_atomic_coordinates"]
