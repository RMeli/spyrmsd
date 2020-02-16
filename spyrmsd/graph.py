import warnings
from typing import Any, Dict, List, Optional, Union

import graph_tool as gt
from graph_tool import topology
import numpy as np
import qcelemental as qcel

# TODO: Move elsewhere?
connectivity_tolerance: float = 0.4


def graph_from_adjacency_matrix(
    adjacency_matrix: Union[np.ndarray, List[List[int]]],
    atomicnums: Optional[Union[np.ndarray, List[int]]] = None,
) -> gt.Graph:
    """
    Graph from adjacency matrix.

    Parameters
    ----------
    adjacency_matrix: Union[np.ndarray, List[List[int]]]
        Adjacency matrix
    atomicnums: Union[np.ndarray, List[int]], optional
        Atomic numbers

    Returns
    -------
    gt.Graph
        graph-tool graph

    Notes
    -----
    It the atomic numbers are passed, they are used as node attributes.
    """

    # Get upper triangular adjacency matrix
    adj = np.triu(adjacency_matrix)

    G = gt.Graph(directed=False)
    G.add_edge_list(np.transpose(adj.nonzero()))

    if atomicnums is not None:
        vprop = G.new_vertex_property("short") # Create property map
        vprop.a = atomicnums # Assign atomic numbers to property map array
        G.vertex_properties["atomicnum"] = vprop # Set property map
    return G


def adjacency_matrix_from_atomic_coordinates(
    atomicnums: np.ndarray, coordinates: np.ndarray
) -> np.ndarray:
    """
    Compute adjacency matrix from atomic coordinates.

    Parameters
    ----------
    atomicnums: numpy.ndarray
        Atomic numbers
    coordinates: numpy.ndarray
        Atomic coordinates

    Returns
    -------
    numpy.ndarray
        Adjacency matrix

    Notes
    -----

    This function is based on an automatic bond perception algorithm: two
    atoms are considered to be bonded when their distance is smaller than
    the sum of their covalent radii plus a tolerance value. [3]_

    .. warning::
        The automatic bond perceptron rule implemented in this functions
        is very simple and only depends on atomic coordinates. Use
        with care!

    .. [3] E. C. Meng and R. A. Lewis, *Determination of molecular topology and atomic
       hybridization states from heavy atom coordinates*, J. Comp. Chem. **12**, 891-898
       (1991).
    """

    n = len(atomicnums)

    assert coordinates.shape == (n, 3)

    A = np.zeros((n, n))

    for i in range(n):
        r_i = qcel.covalentradii.get(atomicnums[i], units="angstrom")

        for j in range(i + 1, n):
            r_j = qcel.covalentradii.get(atomicnums[i], units="angstrom")

            distance = np.sqrt(np.sum((coordinates[i] - coordinates[j]) ** 2))

            if distance < (r_i + r_j + connectivity_tolerance):
                A[i, j] = A[j, i] = 1

    return A


def match_graphs(G1: gt.Graph, G2: gt.Graph) -> List[Dict[Any, Any]]:
    """
    Compute RMSD using the quaternion polynomial method.

    Parameters
    ----------
    G1: networkx.Graph
        Graph 1
    G2: networkx.Graph
        Graph 2

    Returns
    -------
    List[Dict[Any, Any]]
        All possible mappings between nodes of graph 1 and graph 2 (isomorphisms)

    Raises
    ------
    ValueError
        If the graphs `G1` and `G2` are not isomorphic
    """

    try:
        maps = topology.subgraph_isomorphism(G1, G2, vertex_label=(G1.vertex_properties["atomicnum"], G2.vertex_properties["atomicnum"]), subgraph=False)
    except KeyError: # No "atomicnum" vertex property
        warnings.warn(
            "No atomic number information stored on nodes. "
            + "Node matching is not performed..."
        )

        maps = topology.subgraph_isomorphism(G1, G2, subgraph=False)

        print(maps)

    # Check if graphs are actually isomorphic
    if len(maps) == 0:
        # TODO: Create a new exception
        raise ValueError(f"Graphs {G1} and {G2} are not isomorphic.")

    # Extract all isomorphisms in a list
    return [m.a for m in maps]
