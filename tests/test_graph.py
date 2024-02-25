import numpy as np
import pytest

from spyrmsd import constants, graph, io, molecule
from spyrmsd.exceptions import NonIsomorphicGraphs
from spyrmsd.graphs import _common as gc
from spyrmsd.graphs import gt, nx
from tests import molecules

@pytest.mark.parametrize(
    "graph_backend",
    [gt, nx],
)
def test_adjacency_matrix_from_atomic_coordinates_distance(graph_backend) -> None:
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
    G = graph_backend.graph_from_adjacency_matrix(A)

    assert graph_backend.num_edges(G) == 1

@pytest.mark.parametrize(
    "graph_backend",
    [gt, nx],
)
@pytest.mark.parametrize(
    "mol, n_bonds",
    [(molecules.benzene, 12), (molecules.ethanol, 8), (molecules.dialanine, 22)],
)
def test_adjacency_matrix_from_atomic_coordinates(
    mol: molecule.Molecule, n_bonds: int, graph_backend
) -> None:
    A = graph.adjacency_matrix_from_atomic_coordinates(mol.atomicnums, mol.coordinates)

    G = graph_backend.graph_from_adjacency_matrix(A)

    assert graph_backend.num_vertices(G) == len(mol)
    assert graph_backend.num_edges(G) == n_bonds


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

@pytest.mark.parametrize("graph_backend",[gt, nx])
@pytest.mark.parametrize("mol", molecules.allobmolecules)
def test_graph_from_adjacency_matrix(mol, graph_backend) -> None:
    natoms = io.numatoms(mol)
    nbonds = io.numbonds(mol)

    A = io.adjacency_matrix(mol)

    assert A.shape == (natoms, natoms)
    assert np.all(A == A.T)
    assert np.sum(A) == nbonds * 2

    G = graph_backend.graph_from_adjacency_matrix(A)

    assert graph_backend.num_vertices(G) == natoms
    assert graph_backend.num_edges(G) == nbonds

@pytest.mark.parametrize(
    "graph_backend, graph_backend_name",
    [(gt,"gt"), (nx, "nx")],
)
@pytest.mark.parametrize(
    "rawmol, mol", zip(molecules.allobmolecules, molecules.allmolecules)
)
def test_graph_from_adjacency_matrix_atomicnums(rawmol, mol, graph_backend, graph_backend_name) -> None:
    natoms = io.numatoms(rawmol)
    nbonds = io.numbonds(rawmol)

    A = io.adjacency_matrix(rawmol)

    assert len(mol) == natoms
    assert mol.adjacency_matrix.shape == (natoms, natoms)
    assert np.all(mol.adjacency_matrix == A)
    assert np.sum(mol.adjacency_matrix) == nbonds * 2

    G = mol.to_graph(backend=graph_backend_name)

    assert graph_backend.num_vertices(G) == natoms
    assert graph_backend.num_edges(G) == nbonds

    for idx, atomicnum in enumerate(mol.atomicnums):
        assert graph_backend.vertex_property(G, "aprops", idx) == atomicnum


@pytest.mark.parametrize(
    "G1, G2, graph_backend",
    [
        *[(gt.lattice(n, n), gt.lattice(n, n), gt) for n in range(2, 5)],
        *[(nx.lattice(n, n), nx.lattice(n, n), nx) for n in range(2, 5)],
        *[(gt.cycle(n), gt.cycle(n), gt) for n in range(2, 5)],
        *[(nx.cycle(n), nx.cycle(n), nx) for n in range(2, 5)],
    ],
)
def test_match_graphs_isomorphic(G1, G2, graph_backend) -> None:
    with pytest.warns(UserWarning, match=gc.warn_no_atomic_properties):
        isomorphisms = graph_backend.match_graphs(G1, G2)

    assert len(isomorphisms) != 0

@pytest.mark.parametrize(
    "G1, G2, graph_backend",
    [
        *[(gt.lattice(n, n), gt.lattice(n + 1, n), gt) for n in range(2, 5)],
        *[(nx.lattice(n, n), nx.lattice(n + 1, n), nx) for n in range(2, 5)],
        *[(gt.cycle(n), gt.cycle(n + 1), gt) for n in range(1, 5)],
        *[(nx.cycle(n), nx.cycle(n + 1), nx) for n in range(1, 5)],
    ]
)
def test_match_graphs_not_isomorphic(G1, G2, graph_backend) -> None:
    with pytest.raises(
        graph_backend.NonIsomorphicGraphs, match=gc.error_non_isomorphic_graphs
    ), pytest.warns(UserWarning, match=gc.warn_no_atomic_properties):
        graph_backend.match_graphs(G1, G2)

@pytest.mark.parametrize(
    "graph_backend",
    [gt, nx],
)
@pytest.mark.parametrize(
    "property",
    [
        np.array([0, 1, 2], dtype=int),
        np.array([0.1, 1.2, 2.3], dtype=float),
        np.array(["H", "H", "H"], dtype=str),
        np.array(["Csp3", "Csp3", "Csp3"], dtype=str),
        ["LongProperty", "LongPropertyProperty", "LongPropertyProperty"],
    ],
    ids=lambda prop: f"{type(prop).__name__}"
)
def test_build_graph_node_features(graph_backend, property) -> None:
    A = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 1]])
    G = graph_backend.graph_from_adjacency_matrix(A, property)

    assert graph_backend.num_edges(G) == 3


def test_build_graph_node_features_unsupported() -> None:

    A = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 1]])

    property = [True, False, True]

    with pytest.raises(ValueError, match="Unsupported property type:"):
        #Only test graph_tool since NetworkX supports all Python objects as node properties.
        _ = gt.graph_from_adjacency_matrix(A, property)
