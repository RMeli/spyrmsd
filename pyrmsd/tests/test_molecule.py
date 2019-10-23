from pyrmsd import molecule, utils
from pyrmsd.tests import molecules

import qcelemental as qcel

import copy

import numpy as np


def test_load_benzene_xyz():

    mol = molecules.benzene

    assert len(mol.atoms) == 12

    num_h = num_c = 0
    for atom in mol.atoms:
        if atom.atomicnum == 1:
            num_h += 1
        if atom.atomicnum == 6:
            num_c += 1
    assert num_h == num_c == 6


def test_load_ethanol_xyz():

    mol = molecules.ethanol

    assert len(mol.atoms) == 9

    num_h = num_c = num_o = 0
    for atom in mol.atoms:
        if atom.atomicnum == 1:
            num_h += 1
        if atom.atomicnum == 6:
            num_c += 1
        if atom.atomicnum == 8:
            num_o += 1
    assert num_h == 6
    assert num_c == 2
    assert num_o == 1


def test_openbabel_to_molecule_benzene():

    mol = molecules.benzene

    m = molecule.openbabel_to_molecule(mol)

    assert m.atomicnums.shape == (12,)
    assert m.coordinates.shape == (12, 3)


def test_molecule_translate():

    for mol in molecules.xyz:
        mt = molecule.openbabel_to_molecule(mol)

        m = copy.deepcopy(mt)

        t = np.array([0.5, 1.1, -0.1])
        mt.translate(t)

        for tcoord, coord in zip(mt.coordinates, m.coordinates):
            assert np.allclose(tcoord - t, coord)


def test_molecule_rotate_z():

    for mol in molecules.xyz:

        m = molecule.openbabel_to_molecule(mol)

        z_axis = np.array([0, 0, 1])

        for angle in [0, 45, 90]:

            rotated = np.zeros((len(m), 3))
            for i, coord in enumerate(m.coordinates):
                rotated[i] = utils.rotate(coord, angle, z_axis, units="deg")

            m.rotate(angle, z_axis, units="deg")

            assert np.allclose(m.coordinates, rotated)

            # Reset
            m.rotate(-angle, z_axis, units="deg")


def test_molecule_rotate():

    for mol in molecules.xyz:

        m = molecule.openbabel_to_molecule(mol)

        axis = np.random.rand(3)

        for angle in np.random.rand(10) * 180:

            rotated = np.zeros((len(m), 3))
            for i, coord in enumerate(m.coordinates):
                rotated[i] = utils.rotate(coord, angle, axis, units="deg")

            m.rotate(angle, axis, units="deg")

            assert np.allclose(m.coordinates, rotated)

            # Reset
            m.rotate(-angle, axis, units="deg")


def test_molecule_center_of_geometry_benzene():

    mol = molecules.benzene
    m = molecule.openbabel_to_molecule(mol)

    assert np.allclose(m.center_of_geometry(), np.zeros(3))


def test_molecule_center_of_mass_benzene():

    mol = molecules.benzene
    m = molecule.openbabel_to_molecule(mol)

    assert np.allclose(m.center_of_mass(), np.zeros(3))


def test_molecule_center_of_mass_H2():

    atomicnums = [1, 1]
    coordinates = [[0, 0, -1], [0, 0, 1]]

    m = molecule.Molecule(atomicnums, coordinates)

    assert np.allclose(m.center_of_mass(), np.zeros(3))


def test_molecule_center_of_mass_HF():

    atomicnums = [1, 9]
    coordinates = [[0, 0, -1], [0, 0, 1]]

    H_mass = qcel.periodictable.to_mass(1)
    F_mass = qcel.periodictable.to_mass(9)

    z_com = (-H_mass + F_mass) / (H_mass + F_mass)

    m = molecule.Molecule(atomicnums, coordinates)

    assert np.allclose(m.center_of_mass(), np.array([0, 0, z_com]))


def test_molecule_strip_dialanine():

    mol = molecules.dialanine

    m = molecule.openbabel_to_molecule(mol)

    assert len(m) == 23

    m.strip()

    assert len(m) == 11
