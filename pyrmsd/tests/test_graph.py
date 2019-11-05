from pyrmsd import graph, molecule
from pyrmsd.tests import molecules

import numpy as np
import networkx as nx

try:
    from openbabel import openbabel as ob
except ImportError:
    import openbabel as ob

import pytest


@pytest.mark.parametrize(
    "mol, n_bonds",
    [(molecules.benzene, 12), (molecules.ethanol, 8), (molecules.dialanine, 22)],
)
def test_graph_from_atomic_coordinates(mol: molecule.Molecule, n_bonds: int) -> None:

    mol = molecules.dialanine

    G = graph.graph_from_atomic_coordinates(mol.atomicnums, mol.coordinates)

    assert G.number_of_nodes() == len(mol)
    assert G.number_of_edges() == 22


def test_graph_from_atomic_coordinates_named_benzene() -> None:

    mol = molecules.benzene

    G = graph.graph_from_atomic_coordinates(mol.atomicnums, mol.coordinates, named=True)

    assert G.number_of_nodes() == len(mol)
    assert G.number_of_edges() == 12

    nodes = G.nodes()
    for node in nodes:
        element = nodes[node]["element"]

        assert element == "H" or element == "C"


@pytest.mark.parametrize("obmol", molecules.allobmolecules)
def test_adjacency_matrix_from_obmol(obmol) -> None:

    natoms = obmol.OBMol.NumAtoms()
    nbonds = obmol.OBMol.NumBonds()

    A = graph.adjacency_matrix_from_obmol(obmol)

    assert A.shape == (natoms, natoms)
    assert np.alltrue(A == A.T)
    assert np.sum(A) == nbonds * 2

    for bond in ob.OBMolBondIter(obmol.OBMol):
        i = bond.GetBeginAtomIdx() - 1
        j = bond.GetEndAtomIdx() - 1

        assert A[i, j] == 1


@pytest.mark.parametrize("obmol", molecules.allobmolecules)
def test_graph_from_adjacency_matrix(obmol) -> None:

    natoms = obmol.OBMol.NumAtoms()
    nbonds = obmol.OBMol.NumBonds()

    A = graph.adjacency_matrix_from_obmol(obmol)

    assert A.shape == (natoms, natoms)
    assert np.alltrue(A == A.T)
    assert np.sum(A) == nbonds * 2

    G = graph.graph_from_adjacency_matrix(A)

    assert G.number_of_nodes() == natoms
    assert G.number_of_edges() == nbonds


@pytest.mark.parametrize(
    "G1, G2",
    [
        *[(nx.path_graph(n), nx.path_graph(n)) for n in range(5)],
        *[(nx.star_graph(n), nx.star_graph(n)) for n in range(5)],
        *[(nx.cycle_graph(n), nx.cycle_graph(n)) for n in range(1, 5)],
    ],
)
def test_match_graphs_isomorphic(G1: nx.Graph, G2: nx.Graph) -> None:

    isomorphisms = graph.match_graphs(G1, G2)

    assert len(isomorphisms) != 0


@pytest.mark.parametrize(
    "G1, G2",
    [
        *[(nx.path_graph(n), nx.path_graph(n + 1)) for n in range(5)],
        *[(nx.star_graph(n), nx.star_graph(n + 1)) for n in range(5)],
        *[(nx.cycle_graph(n), nx.cycle_graph(n + 1)) for n in range(1, 5)],
    ],
)
def test_match_graphs_not_isomorphic(G1: nx.Graph, G2: nx.Graph) -> None:

    with pytest.raises(ValueError):
        graph.match_graphs(G1, G2)
