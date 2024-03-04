import warnings
from typing import Any, List, Optional, Tuple, Union

import rustworkx as rx
import numpy as np

from spyrmsd.exceptions import NonIsomorphicGraphs
from spyrmsd.graphs._common import (
    error_non_isomorphic_graphs,
    warn_disconnected_graph,
    warn_no_atomic_properties,
)


def graph_from_adjacency_matrix(
    adjacency_matrix: Union[np.ndarray, List[List[int]]],
    aprops: Optional[Union[np.ndarray, List[Any]]] = None,
) -> rx.PyGraph:
    """
    Graph from adjacency matrix.

    Parameters
    ----------
    adjacency_matrix: Union[np.ndarray, List[List[int]]]
        Adjacency matrix
    aprops: Union[np.ndarray, List[Any]], optional
        Atomic properties

    Returns
    -------
    Graph
        Molecular graph

    Notes
    -----
    It the atomic numbers are passed, they are used as node attributes.
    """

    G = rx.PyGraph.from_adjacency_matrix(adjacency_matrix)

    if not rx.is_connected(G):
        warnings.warn(warn_disconnected_graph)

    # if aprops is not None:
    #     attributes = {idx: aprops for idx, aprops in enumerate(aprops)}
    #     nx.set_node_attributes(G, attributes, "aprops")

    return G


def match_graphs(G1, G2) -> List[Tuple[List[int], List[int]]]:
    """
    Compute graph isomorphisms.

    Parameters
    ----------
    G1:
        Graph 1
    G2:
        Graph 2

    Returns
    -------
    List[Tuple[List[int],List[int]]]
        All possible mappings between nodes of graph 1 and graph 2 (isomorphisms)

    Raises
    ------
    NonIsomorphicGraphs
        If the graphs `G1` and `G2` are not isomorphic
    """

    # def match_aprops(node1, node2):
    #     """
    #     Check if atomic properties for two nodes match.
    #     """
    #     return node1["aprops"] == node2["aprops"]

    # if (
    #     nx.get_node_attributes(G1, "aprops") == {}
    #     or nx.get_node_attributes(G2, "aprops") == {}
    # ):
    #     # Nodes without atomic number information
    #     # No node-matching check
    #     node_match = None

    #     warnings.warn(warn_no_atomic_properties)

    # else:
    #     node_match = match_aprops

    #GM = rx.vf2_mapping(G1, G2, node_match)
    GM = rx.vf2_mapping(G1, G2, None)

    # Check if graphs are actually isomorphic
    if not GM:
        raise NonIsomorphicGraphs(error_non_isomorphic_graphs)

    return [
        (list(isomorphism.keys()), list(isomorphism.values()))
        for isomorphism in GM
    ]


def vertex_property(G, vproperty: str, idx: int) -> Any:
    """
    Get vertex (node) property from graph

    Parameters
    ----------
    G:
        Graph
    vproperty: str
        Vertex property name
    idx: int
        Vertex index

    Returns
    -------
    Any
        Vertex property value
    """
    return G.nodes[idx][vproperty]


def num_vertices(G) -> int:
    """
    Number of vertices

    Parameters
    ----------
    G:
        Graph

    Returns
    -------
    int
        Number of vertices (nodes)
    """
    return G.num_nodes()


def num_edges(G) -> int:
    """
    Number of edges

    Parameters
    ----------
    G:
        Graph

    Returns
    -------
    int
        Number of edges
    """
    return G.num_edges()


def lattice(n1, n2):
    """
    Build 2D lattice graph

    Parameters
    ----------
    n1: int
        Number of nodes in dimension 1
    n2: int
        Number of nodes in dimension 2

    Returns
    -------
    Graph
        Lattice graph
    """
    return rx.generators.grid_graph(rows=n1, cols=n2, multigraph=False)


def cycle(n):
    """
    Build cycle graph

    Parameters
    ----------
    n: int
        Number of nodes

    Returns
    -------
    Graph
        Cycle graph
    """
    return rx.generators.cycle_graph(n)
