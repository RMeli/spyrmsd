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


# Results obtained with PyTraj
#   pytraj.analysis.rmsd.rmsd(i, ref=j, nofit=True)
@pytest.mark.parametrize(
    "i, j, result", [(1, 2, 2.60065218), (1, 3, 9.94411523), (2, 3, 9.4091711)]
)
def test_rmsd_dummy_2viz(i: int, j: int, result: float):

    moli = molecule.openbabel_to_molecule(molecules.docking_2viz[i])
    molj = molecule.openbabel_to_molecule(molecules.docking_2viz[j])

    assert rmsd.rmsd_dummy(moli, molj) == pytest.approx(result)


# Results obtained with PyTraj
#   pytraj.analysis.rmsd.rmsd(i, mask="!@H=", ref=j, ref_mask="!@H=", nofit=True)
@pytest.mark.parametrize(
    "i, j, result", [(1, 2, 2.65327362), (1, 3, 10.11099065), (2, 3, 9.57099612)]
)
def test_rmsd_dummy_2viz_stripped(i: int, j: int, result: float):

    moli = molecule.openbabel_to_molecule(molecules.docking_2viz[i])
    molj = molecule.openbabel_to_molecule(molecules.docking_2viz[j])

    moli.strip()
    molj.strip()

    assert rmsd.rmsd_dummy(moli, molj) == pytest.approx(result)


def test_rmsd_dummy_centred_benzene():

    mol1 = molecules.benzene
    mol2 = molecules.benzene

    m1 = molecule.openbabel_to_molecule(mol1)
    m2 = molecule.openbabel_to_molecule(mol2)

    m2.translate(np.array([0, 0, 1]))

    assert rmsd.rmsd_dummy(m1, m2) == pytest.approx(1)
    assert rmsd.rmsd_dummy(m1, m2, center=True) == pytest.approx(0)


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


def test_rmsd_qcp_dialanine():

    mol1 = molecules.dialanine
    mol2 = molecules.dialanine

    m1 = molecule.openbabel_to_molecule(mol1)
    m2 = molecule.openbabel_to_molecule(mol2)

    assert rmsd.rmsd_dummy(m1, m2) == pytest.approx(0)
    assert rmsd.rmsd_qcp(m1, m2) == pytest.approx(0)

    m2.translate(5 * np.random.rand(3))
    m2.rotate(15, np.array([-2, 1, 3]))

    assert rmsd.rmsd_dummy(m1, m2) > 0.0

    assert rmsd.rmsd_qcp(m1, m2) == pytest.approx(0)


# Results obtained with PyTraj
#   pytraj.analysis.rmsd.rmsd(i, ref=j)
@pytest.mark.parametrize(
    "i, j, result", [(1, 2, 1.95277757), (1, 3, 3.11801105), (2, 3, 2.98609758)]
)
def test_rmsd_qcp_2viz(i: int, j: int, result: float):

    moli = molecule.openbabel_to_molecule(molecules.docking_2viz[i])
    molj = molecule.openbabel_to_molecule(molecules.docking_2viz[j])

    assert rmsd.rmsd_qcp(moli, molj) == pytest.approx(result)


# Results obtained with PyTraj
#   pytraj.analysis.rmsd.rmsd(i, "!@H=", ref=j, ref_mask="!@H=")
@pytest.mark.parametrize(
    "i, j, result", [(1, 2, 1.98171656), (1, 3, 3.01799306), (2, 3, 2.82917355)]
)
def test_rmsd_qcp_2viz_stripped(i: int, j: int, result: float):

    moli = molecule.openbabel_to_molecule(molecules.docking_2viz[i])
    molj = molecule.openbabel_to_molecule(molecules.docking_2viz[j])

    # Strip hydrogen atoms
    moli.strip()
    molj.strip()

    assert rmsd.rmsd_qcp(moli, molj) == pytest.approx(result)
