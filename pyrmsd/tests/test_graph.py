from pyrmsd import graph, molecule
from pyrmsd.tests import molecules

import numpy as np

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
