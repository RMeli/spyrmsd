from pyrmsd import rmsd, molecule
from pyrmsd.tests import molecules

import numpy as np

import pytest


def test_rmsd_benzene():

    mol1 = molecules.benzene
    mol2 = molecules.benzene

    m1 = molecule.openbabel_to_molecule(mol1)
    m2 = molecule.openbabel_to_molecule(mol2)

    assert rmsd.rmsd_dummy(m1, m2) == pytest.approx(0)


def test_rmsd_shifted_benzene():

    mol1 = molecules.benzene
    mol2 = molecules.benzene

    m1 = molecule.openbabel_to_molecule(mol1)
    m2 = molecule.openbabel_to_molecule(mol2)

    m2.translate(np.array([0, 0, 1]))

    assert rmsd.rmsd_dummy(m1, m2) == pytest.approx(1)
