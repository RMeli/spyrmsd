import warnings
from typing import Any, Dict, List, Optional, Union

import networkx as nx
import numpy as np

import collections


def graph_from_adjacency_matrix(
    adjacency_matrix: Union[np.ndarray, List[List[int]]],
    atomicnums: Optional[Union[np.ndarray, List[int]]] = None,
) -> nx.Graph:
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
    nx.Graph
        NetworkX graph

    Notes
    -----
    It the atomic numbers are passed, they are used as node attributes.
    """

    G = nx.Graph(adjacency_matrix)

    if atomicnums is not None:
        attributes = {idx: atomicnum for idx, atomicnum in enumerate(atomicnums)}
        nx.set_node_attributes(G, attributes, "atomicnum")

    return G


def match_graphs(G1, G2) -> List[Dict[Any, Any]]:
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

    def match_atomicnum(node1, node2):
        return node1["atomicnum"] == node2["atomicnum"]

    if (
        nx.get_node_attributes(G1, "atomicnum") == {}
        or nx.get_node_attributes(G2, "atomicnum") == {}
    ):
        # Nodes without atomic number information
        # No node-matching check
        node_match = None

        warnings.warn(
            "No atomic number information stored on nodes. "
            + "Node matching is not performed..."
        )

    else:
        node_match = match_atomicnum

    GM = nx.algorithms.isomorphism.GraphMatcher(G1, G2, node_match)

    # Check if graphs are actually isomorphic
    if not GM.is_isomorphic():
        # TODO: Create a new exception
        raise ValueError(f"Graphs {G1} and {G2} are not isomorphic.")

    # Extract all isomorphisms in a list of ordered dicts
    # OrderedDict is necessary to conform with the graph_tool implementation
    # TODO: Improve performance
    isomorphisms = [
        collections.OrderedDict(sorted(isomorphism.items()))
        for isomorphism in GM.isomorphisms_iter()
    ]

    return [np.array(list(isomorphism.values())) for isomorphism in isomorphisms]


def vertex_property(G, vproperty: str, idx: int) -> Any:
    return G.nodes[idx][vproperty]


def num_vertices(G) -> int:
    return G.number_of_nodes()


def num_edges(G) -> int:
    return G.number_of_edges()


def lattice(n1, n2):
    return nx.generators.lattice.grid_2d_graph(n1, n2)


def cycle(n):
    return nx.cycle_graph(n)
