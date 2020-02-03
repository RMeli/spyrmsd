import copy
import os
from collections import defaultdict
from typing import DefaultDict, List, Tuple

import numpy as np
import pytest
import qcelemental as qcel

from spyrmsd import io, molecule, utils
from tests import molecules


# atoms is a list of atomic numbers and atom counts
@pytest.mark.parametrize(
    "mol, atoms",
    [
        (molecules.benzene, [(1, 6), (6, 6)]),
        (molecules.ethanol, [(1, 6), (6, 2), (8, 1)]),
        (molecules.dialanine, [(1, 12), (6, 6), (7, 2), (8, 3)]),
    ],
)
def test_load(mol: molecule.Molecule, atoms: List[Tuple[int, int]]) -> None:

    n = sum([n_atoms for _, n_atoms in atoms])

    assert len(mol) == n
    assert mol.atomicnums.shape == (n,)
    assert mol.coordinates.shape == (n, 3)

    # Count number of atoms of different elements
    atomcount: DefaultDict[int, int] = defaultdict(int)
    for atomicnum in mol.atomicnums:
        atomcount[atomicnum] += 1

    assert len(atomcount) == len(atoms)

    for Z, n_atoms in atoms:
        assert atomcount[Z] == n_atoms


def test_loadall() -> None:

    path = os.path.join(molecules.molpath, "1cbr_docking.sdf")

    mols = io.loadall(path)

    assert len(mols) == 10


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_molecule_translate(mol: molecule.Molecule) -> None:

    mt = copy.deepcopy(mol)

    t = np.array([0.5, 1.1, -0.1])
    mt.translate(t)

    for tcoord, coord in zip(mt.coordinates, mol.coordinates):
        assert np.allclose(tcoord - t, coord)


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_molecule_rotate_z(mol: molecule.Molecule) -> None:

    z_axis = np.array([0, 0, 1])

    for angle in [0, 45, 90]:

        rotated = np.zeros((len(mol), 3))
        for i, coord in enumerate(mol.coordinates):
            rotated[i] = utils.rotate(coord, angle, z_axis, units="deg")

        mol.rotate(angle, z_axis, units="deg")

        assert np.allclose(mol.coordinates, rotated)

        # Reset
        mol.rotate(-angle, z_axis, units="deg")


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_molecule_rotate(mol: molecule.Molecule) -> None:

    axis = np.random.rand(3)

    for angle in np.random.rand(10) * 180:

        rotated = np.zeros((len(mol), 3))
        for i, coord in enumerate(mol.coordinates):
            rotated[i] = utils.rotate(coord, angle, axis, units="deg")

        mol.rotate(angle, axis, units="deg")

        assert np.allclose(mol.coordinates, rotated)

        # Reset
        mol.rotate(-angle, axis, units="deg")


def test_molecule_center_of_geometry_benzene() -> None:

    mol = molecules.benzene

    assert np.allclose(mol.center_of_geometry(), np.zeros(3))


def test_molecule_center_of_mass_benzene() -> None:

    mol = molecules.benzene

    assert np.allclose(mol.center_of_mass(), np.zeros(3))


def test_molecule_center_of_mass_H2() -> None:

    atomicnums = [1, 1]
    coordinates = [[0, 0, -1], [0, 0, 1]]

    mol = molecule.Molecule(atomicnums, coordinates)

    assert np.allclose(mol.center_of_mass(), np.zeros(3))


def test_molecule_center_of_mass_HF() -> None:

    atomicnums = [1, 9]
    coordinates = [[0, 0, -1], [0, 0, 1]]

    H_mass = qcel.periodictable.to_mass(1)
    F_mass = qcel.periodictable.to_mass(9)

    z_com = (-H_mass + F_mass) / (H_mass + F_mass)

    mol = molecule.Molecule(atomicnums, coordinates)

    assert np.allclose(mol.center_of_mass(), np.array([0, 0, z_com]))


@pytest.mark.parametrize(
    "mol, n_atoms, stripped",
    [
        (molecules.benzene, 12, 6),
        (molecules.ethanol, 9, 6),
        (molecules.dialanine, 23, 12),
    ],
)
def test_molecule_strip(mol: molecule.Molecule, n_atoms: int, stripped: int) -> None:

    m = copy.deepcopy(mol)

    assert len(m) == n_atoms

    m.strip()

    assert len(m) == n_atoms - stripped


@pytest.mark.parametrize(
    "mol, n_bonds",
    [(molecules.benzene, 12), (molecules.ethanol, 8), (molecules.dialanine, 22)],
)
def test_graph_from_adjacency_matrix(mol: molecule.Molecule, n_bonds: int) -> None:

    G = mol.to_graph()

    assert G.number_of_nodes() == len(mol)
    assert G.number_of_edges() == n_bonds

    for idx, atomicnum in enumerate(mol.atomicnums):
        assert G.nodes[idx]["atomicnum"] == atomicnum


@pytest.mark.parametrize(
    "mol, n_bonds",
    [(molecules.benzene, 12), (molecules.ethanol, 8), (molecules.dialanine, 22)],
)
def test_graph_from_atomic_coordinates_perception(
    mol: molecule.Molecule, n_bonds: int
) -> None:

    m = copy.deepcopy(mol)

    m.adjacency_matrix = None
    m.G = None

    with pytest.warns(UserWarning):

        # Uses automatic bond perception
        G = m.to_graph()

        assert G.number_of_nodes() == len(m)
        assert G.number_of_edges() == n_bonds

        for idx, atomicnum in enumerate(mol.atomicnums):
            assert G.nodes[idx]["atomicnum"] == atomicnum
