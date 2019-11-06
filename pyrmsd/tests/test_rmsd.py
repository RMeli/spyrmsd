from pyrmsd import rmsd, molecule
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


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_rmsd_qcp(mol: molecule.Molecule) -> None:

    mol1 = copy.deepcopy(mol)
    mol2 = copy.deepcopy(mol)

    assert rmsd.rmsd_dummy(mol1, mol2) == pytest.approx(0)
    assert rmsd.rmsd_qcp(mol1, mol2) == pytest.approx(0)

    for _ in range(10):
        mol2.translate(np.random.rand(3))
        mol2.rotate(np.random.rand(1), np.random.rand(3))

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


@pytest.mark.parametrize(
    "angle, tol", [(60, 1e-5), (120, 1e-5), (180, 1e-12), (240, 1e-5), (300, 1e-5)]
)
def test_rmsd_hungarian_benzene_rotated(angle: float, tol: float) -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    assert rmsd.rmsd_dummy(mol1, mol2) == pytest.approx(0)
    assert rmsd.rmsd_hungarian(mol1, mol2) == pytest.approx(0)

    # Rotations different than 180 degrees introduce numerical errors (~1e-6)
    mol2.rotate(angle, [0, 0, 1], units="deg")

    assert rmsd.rmsd_dummy(mol1, mol2) > 0
    assert rmsd.rmsd_hungarian(mol1, mol2) == pytest.approx(0, abs=tol)


@pytest.mark.parametrize(
    "angle, tol", [(60, 1e-10), (120, 1e-9), (180, 1e-12), (240, 1e-9), (300, 1e-9)]
)
def test_rmsd_hungarian_benzene_shifted_rotated(angle: float, tol: float) -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.translate([0, 0, 1])

    assert rmsd.rmsd_dummy(mol1, mol2) == pytest.approx(1)
    assert rmsd.rmsd_hungarian(mol1, mol2) == pytest.approx(1)

    # Rotations different than 180 degrees introduce numerical errors (~1e-11)
    mol2.rotate(angle, [0, 0, 1], units="deg")

    assert rmsd.rmsd_dummy(mol1, mol2) > 1
    assert rmsd.rmsd_hungarian(mol1, mol2) == pytest.approx(1, abs=tol)


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_rmsd_hungarian_centred(mol: molecule.Molecule) -> None:

    mol1 = copy.deepcopy(mol)
    mol2 = copy.deepcopy(mol)

    mol2.translate(np.random.rand(3))

    assert rmsd.rmsd_hungarian(mol1, mol2) > 0
    assert rmsd.rmsd_hungarian(mol1, mol2, center=True) == pytest.approx(0)


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_rmsd_isomorphic_centred(mol: molecule.Molecule) -> None:

    mol1 = copy.deepcopy(mol)
    mol2 = copy.deepcopy(mol)

    mol2.translate(np.random.rand(3))

    assert rmsd.rmsd_isomorphic(mol1, mol2) > 0
    assert rmsd.rmsd_isomorphic(mol1, mol2, center=True) == pytest.approx(0)


@pytest.mark.parametrize("angle", [60, 120, 180, 240, 300, 360])
def test_rmsd_isomorphic_rotated_benzene(angle: float) -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.rotate(angle, np.array([0, 0, 1]), units="deg")

    assert rmsd.rmsd_dummy(mol1, mol2) > 0
    assert rmsd.rmsd_qcp(mol1, mol2) == pytest.approx(0)
    assert rmsd.rmsd_hungarian(mol1, mol2) == pytest.approx(0, abs=1e-5)
    assert rmsd.rmsd_isomorphic(mol1, mol2) == pytest.approx(0, abs=1e-5)


@pytest.mark.parametrize("angle", [60, 120, 180, 240, 300, 360])
def test_rmsd_isomorphic_rotated_benzene_stripped(angle: float) -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.rotate(angle, np.array([0, 0, 1]), units="deg")

    mol1.strip()
    mol2.strip()

    assert rmsd.rmsd_dummy(mol1, mol2) > 0
    assert rmsd.rmsd_qcp(mol1, mol2) == pytest.approx(0)
    assert rmsd.rmsd_hungarian(mol1, mol2) == pytest.approx(0, abs=1e-5)
    assert rmsd.rmsd_isomorphic(mol1, mol2) == pytest.approx(0, abs=1e-5)
