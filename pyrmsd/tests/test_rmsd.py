from pyrmsd import rmsd
from pyrmsd.tests import molecules

import numpy as np
import copy

import pytest


def test_rmsd_dummy_benzene() -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    assert rmsd.rmsd_dummy(mol1, mol2) == pytest.approx(0)


def test_rmsd_dummy_shifted_benzene() -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.translate(np.array([0, 0, 1]))

    assert rmsd.rmsd_dummy(mol1, mol2) == pytest.approx(1)


# Results obtained with PyTraj
#   pytraj.analysis.rmsd.rmsd(i, ref=j, nofit=True)
@pytest.mark.parametrize(
    "i, j, result", [(1, 2, 2.60065218), (1, 3, 9.94411523), (2, 3, 9.4091711)]
)
def test_rmsd_dummy_2viz(i: int, j: int, result: float) -> None:

    moli = copy.deepcopy(molecules.docking_2viz[i])
    molj = copy.deepcopy(molecules.docking_2viz[j])

    assert rmsd.rmsd_dummy(moli, molj) == pytest.approx(result)


# Results obtained with PyTraj
#   pytraj.analysis.rmsd.rmsd(i, mask="!@H=", ref=j, ref_mask="!@H=", nofit=True)
@pytest.mark.parametrize(
    "i, j, result", [(1, 2, 2.65327362), (1, 3, 10.11099065), (2, 3, 9.57099612)]
)
def test_rmsd_dummy_2viz_stripped(i: int, j: int, result: float) -> None:

    moli = copy.deepcopy(molecules.docking_2viz[i])
    molj = copy.deepcopy(molecules.docking_2viz[j])

    moli.strip()
    molj.strip()

    assert rmsd.rmsd_dummy(moli, molj) == pytest.approx(result)


def test_rmsd_dummy_centred_benzene() -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.translate(np.array([0, 0, 1]))

    assert rmsd.rmsd_dummy(mol1, mol2) == pytest.approx(1)
    assert rmsd.rmsd_dummy(mol1, mol2, center=True) == pytest.approx(0)


def test_rmsd_qcp_benzene() -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    assert rmsd.rmsd_dummy(mol1, mol2) == pytest.approx(0)
    assert rmsd.rmsd_qcp(mol1, mol2) == pytest.approx(0)

    mol2.rotate(10, np.array([1, -2, -1]))

    assert rmsd.rmsd_dummy(mol1, mol2) > 0.0

    assert rmsd.rmsd_qcp(mol1, mol2) == pytest.approx(0)


def test_rmsd_qcp_ethanol() -> None:

    mol1 = copy.deepcopy(molecules.ethanol)
    mol2 = copy.deepcopy(molecules.ethanol)

    assert rmsd.rmsd_dummy(mol1, mol2) == pytest.approx(0)
    assert rmsd.rmsd_qcp(mol1, mol2) == pytest.approx(0)

    cog = mol2.center_of_geometry()
    mol2.translate(-cog)
    mol2.rotate(15, np.array([-2, 1, 3]))
    mol2.translate(cog)

    assert rmsd.rmsd_dummy(mol1, mol2) > 0.0

    assert rmsd.rmsd_qcp(mol1, mol2) == pytest.approx(0)


def test_rmsd_qcp_dialanine() -> None:

    mol1 = copy.deepcopy(molecules.dialanine)
    mol2 = copy.deepcopy(molecules.dialanine)

    assert rmsd.rmsd_dummy(mol1, mol2) == pytest.approx(0)
    assert rmsd.rmsd_qcp(mol1, mol2) == pytest.approx(0)

    mol2.translate(5 * np.random.rand(3))
    mol2.rotate(15, np.array([-2, 1, 3]))

    assert rmsd.rmsd_dummy(mol1, mol2) > 0.0

    assert rmsd.rmsd_qcp(mol1, mol2) == pytest.approx(0)


# Results obtained with PyTraj
#   pytraj.analysis.rmsd.rmsd(i, ref=j)
@pytest.mark.parametrize(
    "i, j, result", [(1, 2, 1.95277757), (1, 3, 3.11801105), (2, 3, 2.98609758)]
)
def test_rmsd_qcp_2viz(i: int, j: int, result: float) -> None:

    moli = copy.deepcopy(molecules.docking_2viz[i])
    molj = copy.deepcopy(molecules.docking_2viz[j])

    assert rmsd.rmsd_qcp(moli, molj) == pytest.approx(result)


# Results obtained with PyTraj
#   pytraj.analysis.rmsd.rmsd(i, "!@H=", ref=j, ref_mask="!@H=")
@pytest.mark.parametrize(
    "i, j, result", [(1, 2, 1.98171656), (1, 3, 3.01799306), (2, 3, 2.82917355)]
)
def test_rmsd_qcp_2viz_stripped(i: int, j: int, result: float) -> None:

    moli = copy.deepcopy(molecules.docking_2viz[i])
    molj = copy.deepcopy(molecules.docking_2viz[j])

    # Strip hydrogen atoms
    moli.strip()
    molj.strip()

    assert rmsd.rmsd_qcp(moli, molj) == pytest.approx(result)
