import copy
import os
from collections import defaultdict
from typing import DefaultDict

import numpy as np
import pytest

import spyrmsd
from spyrmsd import constants, graph, io, molecule, utils


def test_load(mol) -> None:
    # Atoms for each type
    atoms = {
        "benzene": [(1, 6), (6, 6)],
        "ethanol": [(1, 6), (6, 2), (8, 1)],
        "pyridine": [(1, 5), (6, 5), (7, 1)],
        "dialanine": [(1, 12), (6, 6), (7, 2), (8, 3)],
    }

    n = sum([n_atoms for _, n_atoms in atoms[mol.name]])

    assert mol.n_atoms == n
    assert mol.mol.atomicnums.shape == (n,)
    assert mol.mol.coordinates.shape == (n, 3)

    # Count number of atoms of different elements
    atomcount: DefaultDict[int, int] = defaultdict(int)
    for atomicnum in mol.mol.atomicnums:
        atomcount[atomicnum] += 1

    assert len(atomcount) == len(atoms[mol.name])

    for Z, n_atoms in atoms[mol.name]:
        assert atomcount[Z] == n_atoms


def test_loadall(molpath) -> None:
    path = os.path.join(molpath, "1cbr_docking.sdf")

    mols = io.loadall(path)

    assert len(mols) == 10


def test_molecule_translate(mol) -> None:
    mt = copy.deepcopy(mol.mol)

    t = np.array([0.5, 1.1, -0.1])
    mt.translate(t)

    for tcoord, coord in zip(mt.coordinates, mol.mol.coordinates):
        assert np.allclose(tcoord - t, coord)


def test_molecule_rotate_z(mol) -> None:
    z_axis = np.array([0, 0, 1])

    for angle in [0, 45, 90]:
        rotated = np.zeros((mol.n_atoms, 3))
        for i, coord in enumerate(mol.mol.coordinates):
            rotated[i] = utils.rotate(coord, angle, z_axis, units="deg")

        mol.mol.rotate(angle, z_axis, units="deg")

        assert np.allclose(mol.mol.coordinates, rotated)


def test_molecule_rotate(mol) -> None:
    axis = np.random.rand(3)

    for angle in np.random.rand(10) * 180:
        rotated = np.zeros((mol.n_atoms, 3))
        for i, coord in enumerate(mol.mol.coordinates):
            rotated[i] = utils.rotate(coord, angle, axis, units="deg")

        mol.mol.rotate(angle, axis, units="deg")

        assert np.allclose(mol.mol.coordinates, rotated)


def test_molecule_center_of_geometry_benzene(benzene) -> None:
    assert np.allclose(benzene.mol.center_of_geometry(), np.zeros(3))


def test_molecule_center_of_mass_benzene(benzene) -> None:
    assert np.allclose(benzene.mol.center_of_mass(), np.zeros(3))


def test_molecule_center_of_mass_H2() -> None:
    atomicnums = [1, 1]
    coordinates = [[0.0, 0.0, -1.0], [0.0, 0.0, 1.0]]

    mol = molecule.Molecule(atomicnums, coordinates)

    assert np.allclose(mol.center_of_mass(), np.zeros(3))


def test_molecule_center_of_mass_HF() -> None:
    atomicnums = [1, 9]
    coordinates = [[0.0, 0.0, -1.0], [0.0, 0.0, 1.0]]

    H_mass = constants.anum_to_mass[1]
    F_mass = constants.anum_to_mass[9]

    z_com = (-H_mass + F_mass) / (H_mass + F_mass)

    mol = molecule.Molecule(atomicnums, coordinates)

    assert np.allclose(mol.center_of_mass(), np.array([0, 0, z_com]))


def test_molecule_strip(mol) -> None:
    m = copy.deepcopy(mol.mol)

    assert len(m) == mol.n_atoms

    m.strip()

    assert len(m) == mol.n_atoms - mol.n_h


def test_graph_from_adjacency_matrix(mol) -> None:
    G = mol.mol.to_graph()

    assert graph.num_vertices(G) == mol.n_atoms
    assert graph.num_edges(G) == mol.n_bonds

    for idx, atomicnum in enumerate(mol.mol.atomicnums):
        assert graph.vertex_property(G, "aprops", idx) == atomicnum


def test_graph_from_atomic_coordinates_perception(mol) -> None:
    m = copy.deepcopy(mol.mol)

    delattr(m, "adjacency_matrix")
    m.G = {}

    with pytest.warns(UserWarning):
        # Uses automatic bond perception
        G = m.to_graph()

        assert graph.num_vertices(G) == mol.n_atoms
        assert graph.num_edges(G) == mol.n_bonds

        for idx, atomicnum in enumerate(mol.mol.atomicnums):
            assert graph.vertex_property(G, "aprops", idx) == atomicnum


@pytest.mark.parametrize(
    "adjacency",
    [True, False],
    ids=["adjacency", "no_adjacency"],
)
def test_from_obmol(molpath, adjacency):
    pytest.importorskip("openbabel")

    from spyrmsd.optional import obabel as ob

    # Load molecules with OpenBabel
    path = os.path.join(molpath, "1cbr_docking.sdf")
    mols = ob.loadall(path)

    # Convert OpenBabel molecules to spyrmsd molecules
    mols = [molecule.Molecule.from_obabel(mol, adjacency) for mol in mols]

    assert len(mols) == 10

    for mol in mols:
        assert isinstance(mol, molecule.Molecule)

        if adjacency:
            assert mol.adjacency_matrix is not None
        else:
            with pytest.raises(AttributeError):
                # No adjacency_matrix attribute
                mol.adjacency_matrix


@pytest.mark.parametrize(
    "adjacency",
    [True, False],
    ids=["adjacency", "no_adjacency"],
)
def test_from_rdmol(molpath, adjacency):
    pytest.importorskip("rdkit")

    from spyrmsd.optional import rdkit as rd

    # Load molecules with RDKit
    path = os.path.join(molpath, "1cbr_docking.sdf")
    mols = rd.loadall(path)

    # Convert OpenBabel molecules to spyrmsd molecules
    mols = [molecule.Molecule.from_rdkit(mol, adjacency) for mol in mols]

    assert len(mols) == 10

    for mol in mols:
        assert isinstance(mol, molecule.Molecule)

        if adjacency:
            assert mol.adjacency_matrix is not None
        else:
            with pytest.raises(AttributeError):
                # No adjacency_matrix attribute
                mol.adjacency_matrix


@pytest.mark.skipif(
    # Run test if all supported backends are installed
    not set(spyrmsd.graph._supported_backends) <= set(spyrmsd.available_backends),
    reason="Not all of the required backends are installed",
)
def test_molecule_graph_cache(mol) -> None:
    import graph_tool as gt
    import networkx as nx
    import rustworkx as rx

    ## Graph cache persists from previous tests, manually reset them
    mol.G = {}
    spyrmsd.set_backend("networkx")
    mol.to_graph()

    assert isinstance(mol.G["networkx"], nx.Graph)
    assert "graph_tool" not in mol.G.keys()

    spyrmsd.set_backend("graph-tool")
    mol.to_graph()

    spyrmsd.set_backend("rustworkx")
    mol.to_graph()

    ## Make sure all backends (still) have a cache
    assert isinstance(mol.G["networkx"], nx.Graph)
    assert isinstance(mol.G["graph_tool"], gt.Graph)
    assert isinstance(mol.G["rustworkx"], rx.PyGraph)

    ## Strip the molecule to ensure the cache is reset
    mol.strip()

    assert len(mol.G.items()) == 0
