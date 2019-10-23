from pyrmsd import molecule, utils
from pyrmsd.tests import molecules

import qcelemental as qcel

import copy

import numpy as np


def test_load_benzene() -> None:

    mol = molecules.benzene

    assert len(mol) == 12

    num_h = num_c = 0
    for atomicnum in mol.atomicnums:
        if atomicnum == 1:
            num_h += 1
        if atomicnum == 6:
            num_c += 1
    assert num_h == num_c == 6

    assert mol.atomicnums.shape == (12,)
    assert mol.coordinates.shape == (12, 3)


def test_load_ethanol() -> None:

    mol = molecules.ethanol

    assert len(mol) == 9

    num_h = num_c = num_o = 0
    for atomicnum in mol.atomicnums:
        if atomicnum == 1:
            num_h += 1
        if atomicnum == 6:
            num_c += 1
        if atomicnum == 8:
            num_o += 1
    assert num_h == 6
    assert num_c == 2
    assert num_o == 1


def test_molecule_translate() -> None:

    for mol in molecules.xyz:

        mt = copy.deepcopy(mol)

        t = np.array([0.5, 1.1, -0.1])
        mt.translate(t)

        for tcoord, coord in zip(mt.coordinates, mol.coordinates):
            assert np.allclose(tcoord - t, coord)


def test_molecule_rotate_z() -> None:

    for mol in molecules.xyz:

        z_axis = np.array([0, 0, 1])

        for angle in [0, 45, 90]:

            rotated = np.zeros((len(mol), 3))
            for i, coord in enumerate(mol.coordinates):
                rotated[i] = utils.rotate(coord, angle, z_axis, units="deg")

            mol.rotate(angle, z_axis, units="deg")

            assert np.allclose(mol.coordinates, rotated)

            # Reset
            mol.rotate(-angle, z_axis, units="deg")


def test_molecule_rotate() -> None:

    for mol in molecules.xyz:

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


def test_molecule_strip_dialanine() -> None:

    mol = molecules.dialanine

    assert len(mol) == 23

    mol.strip()

    assert len(mol) == 11
