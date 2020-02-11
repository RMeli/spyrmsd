import copy
from typing import List

import numpy as np
import pytest

from spyrmsd import molecule, rmsd
from tests import molecules


def test_rmsd_dummy_benzene() -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    assert rmsd.rmsd_standard(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0)


def test_rmsd_dummy_shifted_benzene() -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.translate(np.array([0, 0, 1]))

    assert rmsd.rmsd_standard(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(1)


# Results obtained with PyTraj
#   pytraj.analysis.rmsd.rmsd(i, ref=j, nofit=True)
@pytest.mark.parametrize(
    "i, j, result", [(1, 2, 2.60065218), (1, 3, 9.94411523), (2, 3, 9.4091711)]
)
def test_rmsd_dummy_2viz(i: int, j: int, result: float) -> None:

    moli = copy.deepcopy(molecules.docking_2viz[i])
    molj = copy.deepcopy(molecules.docking_2viz[j])

    assert rmsd.rmsd_standard(
        moli.coordinates, molj.coordinates, moli.atomicnums, molj.atomicnums
    ) == pytest.approx(result)


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

    assert rmsd.rmsd_standard(
        moli.coordinates, molj.coordinates, moli.atomicnums, molj.atomicnums
    ) == pytest.approx(result)


def test_rmsd_dummy_centred_benzene() -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.translate(np.array([0, 0, 1]))

    assert rmsd.rmsd_standard(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(1)
    assert rmsd.rmsd_standard(
        mol1.coordinates,
        mol2.coordinates,
        mol1.atomicnums,
        mol2.atomicnums,
        center=True,
    ) == pytest.approx(0)


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_rmsd_qcp(mol: molecule.Molecule) -> None:

    mol1 = copy.deepcopy(mol)
    mol2 = copy.deepcopy(mol)

    assert rmsd.rmsd_standard(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0)
    assert rmsd.rmsd_qcp(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0)

    for _ in range(10):
        mol2.translate(np.random.rand(3))
        mol2.rotate(np.random.rand(1), np.random.rand(3))

        assert (
            rmsd.rmsd_standard(
                mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
            )
            > 0.0
        )

        assert rmsd.rmsd_qcp(
            mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
        ) == pytest.approx(0)


# Results obtained with PyTraj
#   pytraj.analysis.rmsd.rmsd(i, ref=j)
@pytest.mark.parametrize(
    "i, j, result", [(1, 2, 1.95277757), (1, 3, 3.11801105), (2, 3, 2.98609758)]
)
def test_rmsd_qcp_2viz(i: int, j: int, result: float) -> None:

    moli = copy.deepcopy(molecules.docking_2viz[i])
    molj = copy.deepcopy(molecules.docking_2viz[j])

    assert rmsd.rmsd_qcp(
        moli.coordinates, molj.coordinates, moli.atomicnums, molj.atomicnums
    ) == pytest.approx(result)


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

    assert rmsd.rmsd_qcp(
        moli.coordinates, molj.coordinates, moli.atomicnums, molj.atomicnums
    ) == pytest.approx(result)


# Results obtained with MDAnalysis
#   trp0 = mda.Universe("trp0.pdb")
#   c0 = trp0.coord -trp0.select_atoms("protein").center_of_geometry()
#   for i in [1,2,3,4,5]:
#       trp = mda.Universe(f"trp{i}.pdb")
#       rmsd_dummy = mda.analysis.rms.rmsd(trp.coord, trp0.coord)
#       c = trp.coord - trp.select_atoms("protein").center_of_geometry()
#       _, rmsd_min = align.rotation_matrix(tc, tc0)
#       print(rmsd_dummy, rmsd_min)
@pytest.mark.parametrize(
    "i, rmsd_dummy, rmsd_min",
    [
        (1, 4.812480551076202, 1.6578281551053196),
        (2, 6.772045449820714, 1.7175638492348284),
        (3, 9.344911262612964, 1.5946081072641485),
        (4, 9.772939589989000, 2.1234944939308220),
        (5, 8.901837608843241, 2.4894805175766606),
    ],
)
def test_rmsd_qcp_protein(i: int, rmsd_dummy: float, rmsd_min: float):

    mol0 = copy.deepcopy(molecules.trp[0])
    mol = copy.deepcopy(molecules.trp[i])

    assert rmsd.rmsd_standard(
        mol0.coordinates, mol.coordinates, mol0.atomicnums, mol.atomicnums
    ) == pytest.approx(rmsd_dummy)
    assert rmsd.rmsd_qcp(
        mol0.coordinates, mol.coordinates, mol0.atomicnums, mol.atomicnums
    ) == pytest.approx(rmsd_min)


@pytest.mark.parametrize(
    "angle, tol", [(60, 1e-4), (120, 1e-4), (180, 1e-4), (240, 1e-4), (300, 1e-4)]
)
def test_rmsd_hungarian_benzene_rotated(angle: float, tol: float) -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    assert rmsd.rmsd_standard(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0)
    assert rmsd.rmsd_hungarian(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0)

    # Rotations different than 180 degrees introduce numerical errors (~1e-6)
    mol2.rotate(angle, [0, 0, 1], units="deg")

    assert (
        rmsd.rmsd_standard(
            mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
        )
        > 0
    )
    assert rmsd.rmsd_hungarian(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0, abs=tol)


@pytest.mark.parametrize(
    "angle, tol", [(60, 1e-4), (120, 1e-4), (180, 1e-4), (240, 1e-4), (300, 1e-4)]
)
def test_rmsd_hungarian_benzene_shifted_rotated(angle: float, tol: float) -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.translate([0, 0, 1])

    assert rmsd.rmsd_standard(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(1)
    assert rmsd.rmsd_hungarian(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(1)

    # Rotations different than 180 degrees introduce numerical errors (~1e-11)
    mol2.rotate(angle, [0, 0, 1], units="deg")

    assert (
        rmsd.rmsd_standard(
            mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
        )
        > 1
    )
    assert rmsd.rmsd_hungarian(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(1, abs=tol)


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_rmsd_hungarian_centred(mol: molecule.Molecule) -> None:

    mol1 = copy.deepcopy(mol)
    mol2 = copy.deepcopy(mol)

    mol2.translate(np.random.rand(3))

    assert (
        rmsd.rmsd_hungarian(
            mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
        )
        > 0
    )
    assert rmsd.rmsd_hungarian(
        mol1.coordinates,
        mol2.coordinates,
        mol1.atomicnums,
        mol2.atomicnums,
        center=True,
    ) == pytest.approx(0)


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_rmsd_isomorphic_centred(mol: molecule.Molecule) -> None:

    mol1 = copy.deepcopy(mol)
    mol2 = copy.deepcopy(mol)

    mol2.translate(np.random.rand(3))

    assert (
        rmsd.rmsd_isomorphic(
            mol1.coordinates,
            mol2.coordinates,
            mol1.adjacency_matrix,
            mol2.adjacency_matrix,
            mol1.atomicnums,
            mol2.atomicnums,
        )
        > 0
    )
    assert rmsd.rmsd_isomorphic(
        mol1.coordinates,
        mol2.coordinates,
        mol1.adjacency_matrix,
        mol2.adjacency_matrix,
        mol1.atomicnums,
        mol2.atomicnums,
        center=True,
    ) == pytest.approx(0)


@pytest.mark.parametrize("angle", [60, 120, 180, 240, 300, 360])
def test_rmsd_isomorphic_rotated_benzene(angle: float) -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.rotate(angle, np.array([0, 0, 1]), units="deg")

    assert (
        rmsd.rmsd_standard(
            mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
        )
        > 0
    )
    assert rmsd.rmsd_qcp(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0)
    assert rmsd.rmsd_hungarian(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0, abs=1e-4)
    assert rmsd.rmsd_isomorphic(
        mol1.coordinates,
        mol2.coordinates,
        mol1.adjacency_matrix,
        mol2.adjacency_matrix,
        mol1.atomicnums,
        mol2.atomicnums,
    ) == pytest.approx(0, abs=1e-4)


@pytest.mark.parametrize("angle", [60, 120, 180, 240, 300, 360])
def test_rmsd_isomorphic_rotated_benzene_stripped(angle: float) -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.rotate(angle, np.array([0, 0, 1]), units="deg")

    mol1.strip()
    mol2.strip()

    assert (
        rmsd.rmsd_standard(
            mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
        )
        > 0
    )
    assert rmsd.rmsd_qcp(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0)
    assert rmsd.rmsd_hungarian(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0, abs=1e-4)
    assert rmsd.rmsd_isomorphic(
        mol1.coordinates,
        mol2.coordinates,
        mol1.adjacency_matrix,
        mol2.adjacency_matrix,
        mol1.atomicnums,
        mol2.atomicnums,
    ) == pytest.approx(0, abs=1e-4)


def test_rmsd_isomorphic_atomicnums_matching_pyridine_stripped() -> None:

    mol1 = copy.deepcopy(molecules.pyridine)
    mol2 = copy.deepcopy(molecules.pyridine)

    mol2.rotate(60, np.array([0, 0, 1]), units="deg")

    mol1.strip()
    mol2.strip()

    # Standard RMSD, correct in this case
    RMSD = rmsd.rmsd_standard(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    )

    # Isomorphic RMSD without atomic number matching is wrong
    with pytest.warns(UserWarning):
        assert rmsd.rmsd_isomorphic(
            mol1.coordinates,
            mol2.coordinates,
            mol1.adjacency_matrix,
            mol2.adjacency_matrix,
        ) == pytest.approx(0, abs=1e-4)

    # Isomorphic RMSD with atomic number matching is correct
    assert rmsd.rmsd_isomorphic(
        mol1.coordinates,
        mol2.coordinates,
        mol1.adjacency_matrix,
        mol2.adjacency_matrix,
        mol1.atomicnums,
        mol2.atomicnums,
    ) == pytest.approx(RMSD, abs=1e-4)


# Results obtained with OpenBabel
@pytest.mark.parametrize(
    "index, RMSD, minimize",
    [
        (1, 0.592256, False),
        (2, 2.11545, False),
        (3, 2.29824, False),
        (4, 9.45773, False),
        (5, 1.35005, False),
        (6, 9.44356, False),
        (7, 9.59758, False),
        (8, 9.55076, False),
        (9, 2.44067, False),
        (10, 9.6171, False),
        (1, 0.476858, True),
        (2, 1.68089, True),
        (3, 1.50267, True),
        (4, 1.90623, True),
        (5, 1.01324, True),
        (6, 1.31716, True),
        (7, 1.11312, True),
        (8, 1.06044, True),
        (9, 0.965387, True),
        (10, 1.37842, True),
    ],
)
def test_rmsd_isomorphic(index: int, RMSD: float, minimize: bool) -> None:

    molc = copy.deepcopy(molecules.docking_1cbr[0])
    mol = copy.deepcopy(molecules.docking_1cbr[index])

    molc.strip()
    mol.strip()

    assert rmsd.rmsd_isomorphic(
        molc.coordinates,
        mol.coordinates,
        molc.adjacency_matrix,
        mol.adjacency_matrix,
        molc.atomicnums,
        mol.atomicnums,
        minimize=minimize,
    ) == pytest.approx(RMSD, abs=1e-5)


# Results obtained with OpenBabel
@pytest.mark.parametrize(
    "minimize, referenceRMSDs",
    [
        (
            False,
            [
                0.592256,
                2.11545,
                2.29824,
                9.45773,
                1.35005,
                9.44356,
                9.59758,
                9.55076,
                2.44067,
                9.6171,
            ],
        ),
        (
            True,
            [
                0.476858,
                1.68089,
                1.50267,
                1.90623,
                1.01324,
                1.31716,
                1.11312,
                1.06044,
                0.965387,
                1.37842,
            ],
        ),
    ],
)
def test_multirmsd_isomorphic(minimize: bool, referenceRMSDs: List[float]) -> None:

    molc = copy.deepcopy(molecules.docking_1cbr[0])
    mols = [copy.deepcopy(mol) for mol in molecules.docking_1cbr[1:]]

    molc.strip()

    for mol in mols:
        mol.strip()

    RMSDs = rmsd.multirmsd_isomorphic(
        molc.coordinates,
        [mol.coordinates for mol in mols],
        molc.adjacency_matrix,
        mols[0].adjacency_matrix,
        molc.atomicnums,
        mols[0].atomicnums,
        minimize=minimize,
    )

    for RMSD, referenceRMSD in zip(RMSDs, referenceRMSDs):
        assert RMSD == pytest.approx(referenceRMSD, abs=1e-5)


# Results obtained with OpenBabel
@pytest.mark.parametrize(
    "minimize, referenceRMSDs",
    [
        (
            False,
            [
                0.592256,
                2.11545,
                2.29824,
                9.45773,
                1.35005,
                9.44356,
                9.59758,
                9.55076,
                2.44067,
                9.6171,
            ],
        ),
        (
            True,
            [
                0.476858,
                1.68089,
                1.50267,
                1.90623,
                1.01324,
                1.31716,
                1.11312,
                1.06044,
                0.965387,
                1.37842,
            ],
        ),
    ],
)
def test_multirmsd_isomorphic_cache(
    minimize: bool, referenceRMSDs: List[float]
) -> None:

    molc = copy.deepcopy(molecules.docking_1cbr[0])
    mols = [copy.deepcopy(mol) for mol in molecules.docking_1cbr[1:]]

    molc.strip()

    for mol in mols:
        mol.strip()

    RMSDs = rmsd.multirmsd_isomorphic(
        molc.coordinates,
        [mol.coordinates for mol in mols],
        molc.adjacency_matrix,
        mols[0].adjacency_matrix,
        molc.atomicnums,
        mols[0].atomicnums,
        minimize=minimize,
        cache=False,
    )

    for RMSD, referenceRMSD in zip(RMSDs, referenceRMSDs):
        assert RMSD == pytest.approx(referenceRMSD, abs=1e-5)
