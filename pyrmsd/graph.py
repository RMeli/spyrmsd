import networkx as nx
import qcelemental as qcel
import numpy as np

from typing import List, Dict, Any

# TODO: Move elsewhere?
covalent_bond_multiplier: float = 1.2


def graph_from_adjacency_matrix(adjacency_matrix: np.ndarray) -> nx.Graph:
    """
    Graph from andjacency matrix.

    Parameters
    ----------
    adjacency_matrix: np.ndarray
        Adjacency matrix

    Returns
    -------
    nx.Graph
        NetworkX graph
    """

    return nx.convert_matrix.from_numpy_array(adjacency_matrix)


# TODO: Refactor to take a molecule as input
def graph_from_atomic_coordinates(atomicnums, coordinates, named=False):

    n = len(atomicnums)

    assert coordinates.shape == (n, 3)

    G = nx.Graph()

    for i in range(n):
        r_i = qcel.covalentradii.get(atomicnums[i], units="angstrom")

        for j in range(i + 1, n):
            r_j = qcel.covalentradii.get(atomicnums[i], units="angstrom")

            distance = np.sqrt(np.sum((coordinates[i] - coordinates[j]) ** 2))

            if distance < (r_i + r_j) * covalent_bond_multiplier:
                G.add_edge(i, j)

    assert G.number_of_nodes() == n

    if named:
        mapping = {
            i: {"element": qcel.periodictable.to_symbol(anum)}
            for i, anum in enumerate(atomicnums)
        }
        nx.set_node_attributes(G, mapping)

    return G


def match_graphs(G1: nx.Graph, G2: nx.Graph) -> List[Dict[Any, Any]]:
    """
    Compute RMSD using the quaternion polynomial method

    Parameters
    ----------
    G1: networkx.Graph
        Graph 1
    G2: networkx.Graph
        Graph 2

    Raturns
    -------
    List[Dict[Any, Any]]
        All possible mappings between nodes of graph 1 and graph 2 (isomorphisms)

    Raises
    ------
    ValueError
        If the graphs `G1` and `G2` are not isomorphic
    """

    GM = nx.algorithms.isomorphism.GraphMatcher(G1, G2)

    # Check if graphs are actually isomorphic
    if not GM.is_isomorphic():
        # TODO: Create a new exception
        raise ValueError(f"Graphs {G1} and {G2} are not isomorphic.")

    # Extract all isomorphisms in a list
    return [isomorphism for isomorphism in GM.isomorphisms_iter()]


if __name__ == "__main__":

    from pyrmsd import io

    import argparse as ap
    from matplotlib import pyplot as plt

    parser = ap.ArgumentParser(description="Draw molecular graph.")

    parser.add_argument("input", type=str, help="Input file")
    parser.add_argument("-o", "--output", type=str, default=None, help="Input file")

    args = parser.parse_args()

    obmol = io.load(args.input)
    mol = io.openbabel_to_molecule(obmol)

    G = graph_from_atomic_coordinates(mol.atomicnums, mol.coordinates, named=True)
    labels = nx.get_node_attributes(G, "element")

    nx.draw_kamada_kawai(G, labels=labels)
    if args.output is None:
        plt.plot()
    else:
        plt.savefig(args.output)
