from pyrmsd import molecule, utils
from pyrmsd.tests import molecules

import numpy as np


def test_load_benzene():

    mol = molecules.benzene

    assert len(mol.atoms) == 12

    num_h = num_c = 0
    for atom in mol.atoms:
        if atom.atomicnum == 1:
            num_h += 1
        if atom.atomicnum == 6:
            num_c += 1
    assert num_h == num_c == 6


def test_openbabel_to_molecule_benzene():

    mol = molecules.benzene

    m = molecule.openbabel_to_molecule(mol)

    assert m.atomicnums.shape == (12,)
    assert m.coordinates.shape == (12, 3)


def test_molecule_rotate_z():

    mol = molecules.benzene
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
