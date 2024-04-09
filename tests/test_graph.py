import numpy as np
import pytest

import spyrmsd
from spyrmsd import constants, graph, io
from spyrmsd.exceptions import NonIsomorphicGraphs
from spyrmsd.graphs import _common as gc


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


def test_adjacency_matrix_from_atomic_coordinates(mol) -> None:
    A = graph.adjacency_matrix_from_atomic_coordinates(
        mol.mol.atomicnums, mol.mol.coordinates
    )

    G = graph.graph_from_adjacency_matrix(A)

    assert graph.num_vertices(G) == mol.n_atoms
    assert graph.num_edges(G) == mol.n_bonds


def test_adjacency_matrix_from_mol(rawmol) -> None:
    natoms = io.numatoms(rawmol.rawmol)
    nbonds = io.numbonds(rawmol.rawmol)

    assert natoms == rawmol.n_atoms
    assert nbonds == rawmol.n_bonds

    A = io.adjacency_matrix(rawmol.rawmol)

    assert A.shape == (natoms, natoms)
    assert np.all(A == A.T)
    assert np.sum(A) == nbonds * 2

    for i, j in io.bonds(rawmol.rawmol):
        assert A[i, j] == 1


def test_graph_from_adjacency_matrix(rawmol) -> None:
    natoms = io.numatoms(rawmol.rawmol)
    nbonds = io.numbonds(rawmol.rawmol)

    assert natoms == rawmol.n_atoms
    assert nbonds == rawmol.n_bonds

    A = io.adjacency_matrix(rawmol.rawmol)

    assert A.shape == (natoms, natoms)
    assert np.all(A == A.T)
    assert np.sum(A) == nbonds * 2

    G = graph.graph_from_adjacency_matrix(A)

    assert graph.num_vertices(G) == natoms
    assert graph.num_edges(G) == nbonds


def test_graph_from_adjacency_matrix_atomicnums(rawmol) -> None:
    mol = rawmol.mol

    natoms = io.numatoms(rawmol.rawmol)
    nbonds = io.numbonds(rawmol.rawmol)

    A = io.adjacency_matrix(rawmol.rawmol)

    assert len(mol) == natoms
    assert mol.adjacency_matrix.shape == (natoms, natoms)
    assert np.all(mol.adjacency_matrix == A)
    assert np.sum(mol.adjacency_matrix) == nbonds * 2

    G = mol.to_graph()

    assert graph.num_vertices(G) == natoms
    assert graph.num_edges(G) == nbonds

    for idx, atomicnum in enumerate(mol.atomicnums):
        assert graph.vertex_property(G, "aprops", idx) == atomicnum


@pytest.mark.parametrize("n", list(range(2, 5)))
def test_match_graphs_isomorphic_lattice(n) -> None:
    G1 = graph.lattice(n, n)
    G2 = graph.lattice(n, n)

    with pytest.warns(UserWarning, match=gc.warn_no_atomic_properties):
        isomorphisms = graph.match_graphs(G1, G2)

    assert len(isomorphisms) != 0


@pytest.mark.parametrize("n", list(range(2, 5)))
def test_match_graphs_isomorphic_cycle(n) -> None:
    G1 = graph.cycle(n)
    G2 = graph.cycle(n)

    with pytest.warns(UserWarning, match=gc.warn_no_atomic_properties):
        isomorphisms = graph.match_graphs(G1, G2)

    assert len(isomorphisms) != 0


@pytest.mark.parametrize("n", list(range(2, 5)))
def test_match_graphs_not_isomorphic_lattice(n) -> None:
    G1 = graph.lattice(n, n)
    G2 = graph.lattice(n + 1, n)

    with pytest.raises(
        NonIsomorphicGraphs, match=gc.error_non_isomorphic_graphs
    ), pytest.warns(UserWarning, match=gc.warn_no_atomic_properties):
        graph.match_graphs(G1, G2)


@pytest.mark.parametrize("n", range(2, 5))
def test_match_graphs_not_isomorphic_cycle(n) -> None:
    G1 = graph.cycle(n)
    G2 = graph.cycle(n + 1)

    with pytest.raises(
        NonIsomorphicGraphs, match=gc.error_non_isomorphic_graphs
    ), pytest.warns(UserWarning, match=gc.warn_no_atomic_properties):
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
    if spyrmsd.get_backend() != "graph-tool":
        pytest.skip(
            "NetworkX and RustworkX support all Python objects as node properties."
        )

    A = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 1]])

    property = [True, False, True]

    with pytest.raises(ValueError, match="Unsupported property type:"):
        _ = graph.graph_from_adjacency_matrix(A, property)
