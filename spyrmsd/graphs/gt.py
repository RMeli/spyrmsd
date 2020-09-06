import warnings
from typing import Any, List, Optional, Tuple, Union

import graph_tool as gt
import numpy as np
from graph_tool import generation, topology


def graph_from_adjacency_matrix(
    adjacency_matrix: Union[np.ndarray, List[List[int]]],
    atomicnums: Optional[Union[np.ndarray, List[int]]] = None,
):
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
    Graph
        Molecular graph

    Notes
    -----
    It the atomic numbers are passed, they are used as node attributes.
    """

    # Get upper triangular adjacency matrix
    adj = np.triu(adjacency_matrix)

    G = gt.Graph(directed=False)
    G.add_edge_list(np.transpose(adj.nonzero()))

    if atomicnums is not None:
        vprop = G.new_vertex_property("short")  # Create property map (of C type short)
        vprop.a = atomicnums  # Assign atomic numbers to property map array
        G.vertex_properties["atomicnum"] = vprop  # Set property map

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
    ValueError
        If the graphs `G1` and `G2` are not isomorphic
    """

    try:
        maps = topology.subgraph_isomorphism(
            G1,
            G2,
            vertex_label=(
                G1.vertex_properties["atomicnum"],
                G2.vertex_properties["atomicnum"],
            ),
            subgraph=False,
        )
    except KeyError:  # No "atomicnum" vertex property
        warnings.warn(
            "No atomic number information stored on nodes. "
            + "Node matching is not performed..."
        )

        maps = topology.subgraph_isomorphism(G1, G2, subgraph=False)

    # Check if graphs are actually isomorphic
    if len(maps) == 0:
        # TODO: Create a new exception
        raise ValueError(
            "Graphs are not isomorphic."
            "\nMake sure graphs have the same connectivity."
        )

    n = num_vertices(G1)

    # Extract all isomorphisms in a list
    return [(np.arange(0, n, dtype=int), m.a) for m in maps]


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
    return G.vertex_properties[vproperty][idx]


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
    return G.num_vertices()


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


def lattice(n1: int, n2: int):
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
    return generation.lattice((n1, n2))


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
    return generation.circular_graph(n)
