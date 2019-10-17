from pyrmsd import rmsd, molecule
from pyrmsd.tests import molecules

import numpy as np

import pytest


def test_rmsd_dummy_benzene():

    mol1 = molecules.benzene
    mol2 = molecules.benzene

    m1 = molecule.openbabel_to_molecule(mol1)
    m2 = molecule.openbabel_to_molecule(mol2)

    assert rmsd.rmsd_dummy(m1, m2) == pytest.approx(0)


def test_rmsd_dummy_shifted_benzene():

    mol1 = molecules.benzene
    mol2 = molecules.benzene

    m1 = molecule.openbabel_to_molecule(mol1)
    m2 = molecule.openbabel_to_molecule(mol2)

    m2.translate(np.array([0, 0, 1]))

    assert rmsd.rmsd_dummy(m1, m2) == pytest.approx(1)


def test_rmsd_qcp_benzene():

    mol1 = molecules.benzene
    mol2 = molecules.benzene

    m1 = molecule.openbabel_to_molecule(mol1)
    m2 = molecule.openbabel_to_molecule(mol2)

    assert rmsd.rmsd_dummy(m1, m2) == pytest.approx(0)
    assert rmsd.rmsd_qcp(m1, m2) == pytest.approx(0)

    m2.rotate(10, np.array([1, -2, -1]))

    assert rmsd.rmsd_dummy(m1, m2) > 0.0

    assert rmsd.rmsd_qcp(m1, m2) == pytest.approx(0)


def test_rmsd_qcp_ethanol():

    mol1 = molecules.ethanol
    mol2 = molecules.ethanol

    m1 = molecule.openbabel_to_molecule(mol1)
    m2 = molecule.openbabel_to_molecule(mol2)

    assert rmsd.rmsd_dummy(m1, m2) == pytest.approx(0)
    assert rmsd.rmsd_qcp(m1, m2) == pytest.approx(0)

    cog = m2.center_of_geometry()
    m2.translate(-cog)
    m2.rotate(15, np.array([-2, 1, 3]))
    m2.translate(cog)

    assert rmsd.rmsd_dummy(m1, m2) > 0.0

    assert rmsd.rmsd_qcp(m1, m2) == pytest.approx(0)
