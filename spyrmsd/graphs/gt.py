import warnings
from typing import Any, Dict, List, Optional, Union

import graph_tool as gt
from graph_tool import topology
import numpy as np


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
    graph_tool
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
        vprop = G.new_vertex_property("short")  # Create property map
        vprop.a = atomicnums  # Assign atomic numbers to property map array
        G.vertex_properties["atomicnum"] = vprop  # Set property map

    return G


def match_graphs(G1, G2):
    """
    Compute RMSD using the quaternion polynomial method.

    Parameters
    ----------
    G1: graph_tool.Graph
        Graph 1
    G2: graph_tool.Graph
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

        print(maps)

    # Check if graphs are actually isomorphic
    if len(maps) == 0:
        # TODO: Create a new exception
        raise ValueError(f"Graphs {G1} and {G2} are not isomorphic.")

    # Extract all isomorphisms in a list
    return [m.a for m in maps]


def vertex_property(G, vproperty: str, idx: int) -> Any:
    return G.vertex_properties[vproperty][idx]


def num_vertices(G) -> int:
    return G.num_vertices()


def num_edges(G) -> int:
    return G.num_edges()
