import numpy as np
import pytest

from spyrmsd import constants, graph, io, molecule
from spyrmsd.exceptions import NonIsomorphicGraphs
from tests import molecules


def test_adjacency_matrix_from_atomic_coordinates_distance() -> None:
    # Lithium hydride (LiH)
    # H and Li have very different covalent radii
    atomicnums = np.array([1, 3])

    # Distance is the sum of covalent radii
    d = sum([constants.anum_to_covalentradius[anum] for anum in atomicnums])

    # Distance between two atoms is barely enough to create a bond
    # If the covalent radii are not correct, no bond will be created
    coordinates = np.array(
        [[0, 0, 0], [0, 0, d + constants.connectivity_tolerance - 0.01]]
    )

    A = graph.adjacency_matrix_from_atomic_coordinates(atomicnums, coordinates)
    G = graph.graph_from_adjacency_matrix(A)

    assert graph.num_edges(G) == 1


@pytest.mark.parametrize(
    "mol, n_bonds",
    [(molecules.benzene, 12), (molecules.ethanol, 8), (molecules.dialanine, 22)],
)
def test_adjacency_matrix_from_atomic_coordinates(
    mol: molecule.Molecule, n_bonds: int
) -> None:

    A = graph.adjacency_matrix_from_atomic_coordinates(mol.atomicnums, mol.coordinates)

    G = graph.graph_from_adjacency_matrix(A)

    assert graph.num_vertices(G) == len(mol)
    assert graph.num_edges(G) == n_bonds


@pytest.mark.parametrize("mol", molecules.allobmolecules)
def test_adjacency_matrix_from_mol(mol) -> None:

    natoms = io.numatoms(mol)
    nbonds = io.numbonds(mol)

    A = io.adjacency_matrix(mol)

    assert A.shape == (natoms, natoms)
    assert np.all(A == A.T)
    assert np.sum(A) == nbonds * 2

    for i, j in io.bonds(mol):

        assert A[i, j] == 1


@pytest.mark.parametrize("mol", molecules.allobmolecules)
def test_graph_from_adjacency_matrix(mol) -> None:

    natoms = io.numatoms(mol)
    nbonds = io.numbonds(mol)

    A = io.adjacency_matrix(mol)

    assert A.shape == (natoms, natoms)
    assert np.all(A == A.T)
    assert np.sum(A) == nbonds * 2

    G = graph.graph_from_adjacency_matrix(A)

    assert graph.num_vertices(G) == natoms
    assert graph.num_edges(G) == nbonds


@pytest.mark.parametrize(
    "rawmol, mol", zip(molecules.allobmolecules, molecules.allmolecules)
)
def test_graph_from_adjacency_matrix_atomicnums(rawmol, mol) -> None:

    natoms = io.numatoms(rawmol)
    nbonds = io.numbonds(rawmol)

    A = io.adjacency_matrix(rawmol)

    assert len(mol) == natoms
    assert mol.adjacency_matrix.shape == (natoms, natoms)
    assert np.all(mol.adjacency_matrix == A)
    assert np.sum(mol.adjacency_matrix) == nbonds * 2

    G = mol.to_graph()

    assert graph.num_vertices(G) == natoms
    assert graph.num_edges(G) == nbonds

    for idx, atomicnum in enumerate(mol.atomicnums):
        assert graph.vertex_property(G, "aprops", idx) == atomicnum


@pytest.mark.parametrize(
    "G1, G2",
    [
        *[(graph.lattice(n, n), graph.lattice(n, n)) for n in range(2, 5)],
        *[(graph.cycle(n), graph.cycle(n)) for n in range(2, 5)],
    ],
)
def test_match_graphs_isomorphic(G1, G2) -> None:

    with pytest.warns(
        UserWarning, match="No atomic property information stored on nodes."
    ):
        isomorphisms = graph.match_graphs(G1, G2)

    assert len(isomorphisms) != 0


@pytest.mark.parametrize(
    "G1, G2",
    [
        *[(graph.lattice(n, n), graph.lattice(n + 1, n)) for n in range(2, 5)],
        *[(graph.cycle(n), graph.cycle(n + 1)) for n in range(1, 5)],
    ],
)
def test_match_graphs_not_isomorphic(G1, G2) -> None:

    with pytest.raises(
        NonIsomorphicGraphs, match="Graphs are not isomorphic."
    ), pytest.warns(UserWarning, match="No atomic number information stored on nodes."):
        graph.match_graphs(G1, G2)


@pytest.mark.parametrize(
    "property",
    [
        np.array([0, 1, 2], dtype=int),
        np.array([0.1, 1.2, 2.3], dtype=float),
        np.array(["H", "H", "H"], dtype=str),
        np.array(["Csp3", "Csp3", "Csp3"], dtype=str),
        ["LongProperty", "LongPropertyProperty", "LongPropertyProperty"],
    ],
)
def test_build_graph_node_features(property) -> None:
    A = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 1]])
    G = graph.graph_from_adjacency_matrix(A, property)

    assert graph.num_edges(G) == 3


def test_build_graph_node_features_unsupported() -> None:
    pytest.importorskip(
        "graph_tool", reason="NetworkX supports all Python objects as node properties."
    )

    A = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 1]])

    property = [True, False, True]

    with pytest.raises(ValueError, match="Unsupported property type:"):
        _ = graph.graph_from_adjacency_matrix(A, property)
