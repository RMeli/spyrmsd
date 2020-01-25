import networkx as nx
import numpy as np
import pytest

from spyrmsd import graph, io, molecule
from tests import molecules


@pytest.mark.parametrize(
    "mol, n_bonds",
    [(molecules.benzene, 12), (molecules.ethanol, 8), (molecules.dialanine, 22)],
)
def test_adjacency_matrix_from_atomic_coordinates(
    mol: molecule.Molecule, n_bonds: int
) -> None:

    A = graph.adjacency_matrix_from_atomic_coordinates(mol.atomicnums, mol.coordinates)

    G = graph.graph_from_adjacency_matrix(A)

    assert G.number_of_nodes() == len(mol)
    assert G.number_of_edges() == n_bonds


@pytest.mark.parametrize("mol", molecules.allobmolecules)
def test_adjacency_matrix_from_mol(mol) -> None:

    natoms = io.numatoms(mol)
    nbonds = io.numbonds(mol)

    A = io.adjacency_matrix(mol)

    assert A.shape == (natoms, natoms)
    assert np.alltrue(A == A.T)
    assert np.sum(A) == nbonds * 2

    for i, j in io.bonds(mol):

        assert A[i, j] == 1


@pytest.mark.parametrize("mol", molecules.allobmolecules)
def test_graph_from_adjacency_matrix(mol) -> None:

    natoms = io.numatoms(mol)
    nbonds = io.numbonds(mol)

    A = io.adjacency_matrix(mol)

    assert A.shape == (natoms, natoms)
    assert np.alltrue(A == A.T)
    assert np.sum(A) == nbonds * 2

    G = graph.graph_from_adjacency_matrix(A)

    assert G.number_of_nodes() == natoms
    assert G.number_of_edges() == nbonds


@pytest.mark.parametrize(
    "G1, G2",
    [
        *[(nx.path_graph(n), nx.path_graph(n)) for n in range(3)],
        *[(nx.star_graph(n), nx.star_graph(n)) for n in range(3)],
        *[(nx.cycle_graph(n), nx.cycle_graph(n)) for n in range(1, 3)],
    ],
)
def test_match_graphs_isomorphic(G1: nx.Graph, G2: nx.Graph) -> None:

    isomorphisms = graph.match_graphs(G1, G2)

    assert len(isomorphisms) != 0


@pytest.mark.parametrize(
    "G1, G2",
    [
        *[(nx.path_graph(n), nx.path_graph(n + 1)) for n in range(3)],
        *[(nx.star_graph(n), nx.star_graph(n + 1)) for n in range(3)],
        *[(nx.cycle_graph(n), nx.cycle_graph(n + 1)) for n in range(1, 3)],
    ],
)
def test_match_graphs_not_isomorphic(G1: nx.Graph, G2: nx.Graph) -> None:

    with pytest.raises(ValueError):
        graph.match_graphs(G1, G2)
