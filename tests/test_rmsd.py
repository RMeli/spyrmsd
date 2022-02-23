import copy
from typing import List

import numpy as np
import pytest

from spyrmsd import molecule, qcp, rmsd
from tests import molecules


@pytest.fixture(autouse=True, params=[True, False])
def lambda_max_failure(monkeypatch, request):
    """
    Monkey patch fixture for :code:`lambda_max` function to simulate convergence
    failures.

    Notes
    -----
    The :fun:`lambda_max` function can sometimes fail to converge and raises a
    :code:`RuntimeError` (see GitHub Issue #35 by @kjelljorner). If this occours,
    there is an automatic fallback to the explicit calculation of the highest
    eigenvalue. This is not easy to reproduce but need to be tested.

    Using monkey patching we run all tests with the original :fun:`lambda_max` function
    and again with a patched :fun:`lambda_max` function that always raise the
    :code:`RuntimeError`, so that the fallback is tested instead.

    https://github.com/RMeli/spyrmsd/issues/35
    """
    if request.param:
        # Patch lambda_max to always raise an exception
        # This enforces _lambda_max_eig to be used instead
        def lambda_max_failure(Ga, Gb, c2, c1, c0):
            # Simulate Newton method convergence failure
            raise RuntimeError

        monkeypatch.setattr(qcp, "lambda_max", lambda_max_failure)


@pytest.mark.parametrize("t, RMSD", [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)])
def test_rmsd_benzene(t: float, RMSD: float) -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.translate(np.array([0, 0, t]))

    assert rmsd.rmsd(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(RMSD)


# Results obtained with PyTraj
#   pytraj.analysis.rmsd.rmsd(i, ref=j, nofit=True)
@pytest.mark.parametrize(
    "i, j, result", [(1, 2, 2.60065218), (1, 3, 9.94411523), (2, 3, 9.4091711)]
)
def test_rmsd_2viz(i: int, j: int, result: float) -> None:

    moli = copy.deepcopy(molecules.docking_2viz[i])
    molj = copy.deepcopy(molecules.docking_2viz[j])

    assert rmsd.rmsd(
        moli.coordinates, molj.coordinates, moli.atomicnums, molj.atomicnums
    ) == pytest.approx(result)


# Results obtained with PyTraj
#   pytraj.analysis.rmsd.rmsd(i, mask="!@H=", ref=j, ref_mask="!@H=", nofit=True)
@pytest.mark.parametrize(
    "i, j, result", [(1, 2, 2.65327362), (1, 3, 10.11099065), (2, 3, 9.57099612)]
)
def test_rmsd_2viz_stripped(i: int, j: int, result: float) -> None:

    moli = copy.deepcopy(molecules.docking_2viz[i])
    molj = copy.deepcopy(molecules.docking_2viz[j])

    moli.strip()
    molj.strip()

    assert rmsd.rmsd(
        moli.coordinates, molj.coordinates, moli.atomicnums, molj.atomicnums
    ) == pytest.approx(result)


def test_rmsd_centred_benzene() -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.translate(np.array([0, 0, 1]))

    assert rmsd.rmsd(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(1)

    assert (
        rmsd.rmsd(
            mol1.coordinates,
            mol2.coordinates,
            mol1.atomicnums,
            mol2.atomicnums,
            center=True,
        )
        == pytest.approx(0)
    )


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_rmsd_minimize(mol: molecule.Molecule) -> None:

    mol1 = copy.deepcopy(mol)
    mol2 = copy.deepcopy(mol)

    assert rmsd.rmsd(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0)

    assert (
        rmsd.rmsd(
            mol1.coordinates,
            mol2.coordinates,
            mol1.atomicnums,
            mol2.atomicnums,
            minimize=True,
        )
        == pytest.approx(0)
    )

    for _ in range(10):
        mol2.translate(np.random.rand(3))
        mol2.rotate(np.random.rand(1), np.random.rand(3))

        assert (
            rmsd.rmsd(
                mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
            )
            > 0.0
        )

        assert (
            rmsd.rmsd(
                mol1.coordinates,
                mol2.coordinates,
                mol1.atomicnums,
                mol2.atomicnums,
                minimize=True,
            )
            == pytest.approx(0)
        )


# Results obtained with PyTraj
#   pytraj.analysis.rmsd.rmsd(i, ref=j)
@pytest.mark.parametrize(
    "i, j, result", [(1, 2, 1.95277757), (1, 3, 3.11801105), (2, 3, 2.98609758)]
)
def test_rmsd_qcp_2viz(i: int, j: int, result: float) -> None:

    moli = copy.deepcopy(molecules.docking_2viz[i])
    molj = copy.deepcopy(molecules.docking_2viz[j])

    assert (
        rmsd.rmsd(
            moli.coordinates,
            molj.coordinates,
            moli.atomicnums,
            molj.atomicnums,
            minimize=True,
        )
        == pytest.approx(result)
    )


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

    assert (
        rmsd.rmsd(
            moli.coordinates,
            molj.coordinates,
            moli.atomicnums,
            molj.atomicnums,
            minimize=True,
        )
        == pytest.approx(result)
    )


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

    assert rmsd.rmsd(
        mol0.coordinates, mol.coordinates, mol0.atomicnums, mol.atomicnums
    ) == pytest.approx(rmsd_dummy)

    assert (
        rmsd.rmsd(
            mol0.coordinates,
            mol.coordinates,
            mol0.atomicnums,
            mol.atomicnums,
            minimize=True,
        )
        == pytest.approx(rmsd_min)
    )


@pytest.mark.parametrize(
    "angle, tol", [(60, 1e-4), (120, 1e-4), (180, 1e-4), (240, 1e-4), (300, 1e-4)]
)
def test_rmsd_hungarian_benzene_rotated(angle: float, tol: float) -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    assert rmsd.rmsd(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0)

    assert rmsd.hrmsd(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0)

    # Rotations different than 180 degrees introduce numerical errors (~1e-6)
    mol2.rotate(angle, [0, 0, 1], units="deg")

    assert (
        rmsd.rmsd(mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums)
        > 0
    )
    assert rmsd.hrmsd(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0, abs=tol)


@pytest.mark.parametrize("d", [-0.5, 0.0, 0.5, 1.0, 1.5])
@pytest.mark.parametrize(
    "angle, tol", [(60, 1e-4), (120, 1e-4), (180, 1e-4), (240, 1e-4), (300, 1e-4)]
)
def test_rmsd_hungarian_benzene_shifted_rotated(
    d: float, angle: float, tol: float
) -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.translate([0, 0, d])

    assert rmsd.rmsd(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(abs(d))

    assert rmsd.hrmsd(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(abs(d))

    # Rotations different than 180 degrees introduce numerical errors (~1e-11)
    mol2.rotate(angle, [0, 0, 1], units="deg")

    assert rmsd.rmsd(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) > abs(d)
    assert rmsd.hrmsd(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(abs(d), abs=tol)


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_rmsd_hungarian_centred(mol: molecule.Molecule) -> None:

    mol1 = copy.deepcopy(mol)
    mol2 = copy.deepcopy(mol)

    mol2.translate(np.random.rand(3))

    assert (
        rmsd.hrmsd(mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums)
        > 0
    )

    assert (
        rmsd.hrmsd(
            mol1.coordinates,
            mol2.coordinates,
            mol1.atomicnums,
            mol2.atomicnums,
            center=True,
        )
        == pytest.approx(0)
    )


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_symmrmsd_centred(mol: molecule.Molecule) -> None:

    mol1 = copy.deepcopy(mol)
    mol2 = copy.deepcopy(mol)

    mol2.translate(np.random.rand(3))

    assert (
        rmsd.symmrmsd(
            mol1.coordinates,
            mol2.coordinates,
            mol1.atomicnums,
            mol2.atomicnums,
            mol1.adjacency_matrix,
            mol2.adjacency_matrix,
        )
        > 0
    )

    assert (
        rmsd.symmrmsd(
            mol1.coordinates,
            mol2.coordinates,
            mol1.atomicnums,
            mol2.atomicnums,
            mol1.adjacency_matrix,
            mol2.adjacency_matrix,
            center=True,
        )
        == pytest.approx(0)
    )


@pytest.mark.parametrize("angle", [60, 120, 180, 240, 300, 360])
def test_symmrmsd_rotated_benzene(angle: float) -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.rotate(angle, np.array([0, 0, 1]), units="deg")

    assert (
        rmsd.rmsd(mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums)
        > 0
    )

    assert (
        rmsd.rmsd(
            mol1.coordinates,
            mol2.coordinates,
            mol1.atomicnums,
            mol2.atomicnums,
            minimize=True,
        )
        == pytest.approx(0)
    )

    assert rmsd.hrmsd(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0, abs=1e-4)

    assert (
        rmsd.symmrmsd(
            mol1.coordinates,
            mol2.coordinates,
            mol1.atomicnums,
            mol2.atomicnums,
            mol1.adjacency_matrix,
            mol2.adjacency_matrix,
        )
        == pytest.approx(0, abs=1e-4)
    )


@pytest.mark.parametrize("angle", [60, 120, 180, 240, 300, 360])
def test_symmrmsd_rotated_benzene_stripped(angle: float) -> None:

    mol1 = copy.deepcopy(molecules.benzene)
    mol2 = copy.deepcopy(molecules.benzene)

    mol2.rotate(angle, np.array([0, 0, 1]), units="deg")

    mol1.strip()
    mol2.strip()

    assert (
        rmsd.rmsd(mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums)
        > 0
    )

    assert (
        rmsd.rmsd(
            mol1.coordinates,
            mol2.coordinates,
            mol1.atomicnums,
            mol2.atomicnums,
            minimize=True,
        )
        == pytest.approx(0)
    )

    assert rmsd.hrmsd(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    ) == pytest.approx(0, abs=1e-4)

    assert (
        rmsd.symmrmsd(
            mol1.coordinates,
            mol2.coordinates,
            mol1.atomicnums,
            mol2.atomicnums,
            mol1.adjacency_matrix,
            mol2.adjacency_matrix,
        )
        == pytest.approx(0, abs=1e-4)
    )


def test_symmrmsd_atomicnums_matching_pyridine_stripped() -> None:

    mol1 = copy.deepcopy(molecules.pyridine)
    mol2 = copy.deepcopy(molecules.pyridine)

    mol2.rotate(60, np.array([0, 0, 1]), units="deg")

    mol1.strip()
    mol2.strip()

    # Standard RMSD, correct in this case
    RMSD = rmsd.rmsd(
        mol1.coordinates, mol2.coordinates, mol1.atomicnums, mol2.atomicnums
    )

    # Isomorphic RMSD with atomic number matching is correct
    # Without atomic number matching this would be wrong because of higher symmetry
    assert (
        rmsd.symmrmsd(
            mol1.coordinates,
            mol2.coordinates,
            mol1.atomicnums,
            mol2.atomicnums,
            mol1.adjacency_matrix,
            mol2.adjacency_matrix,
        )
        == pytest.approx(RMSD, abs=1e-4)
    )


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
def test_rmsd_symmrmsd(index: int, RMSD: float, minimize: bool) -> None:

    molc = copy.deepcopy(molecules.docking_1cbr[0])
    mol = copy.deepcopy(molecules.docking_1cbr[index])

    molc.strip()
    mol.strip()

    assert (
        rmsd.symmrmsd(
            molc.coordinates,
            mol.coordinates,
            molc.atomicnums,
            mol.atomicnums,
            molc.adjacency_matrix,
            mol.adjacency_matrix,
            minimize=minimize,
        )
        == pytest.approx(RMSD, abs=1e-5)
    )


def test_rmsd_symmrmsd_disconnected_node() -> None:

    c = np.array([[0.0, 1.0, 2.0], [1.0, 2.0, 3.0], [2.0, 3.0, 4.0]])
    a = np.array([0, 1, 4])

    # Adjacency matrix with disconnected node
    A = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])

    with pytest.warns(UserWarning, match="Disconnected graph detected."):
        assert rmsd.symmrmsd(c, c, a, a, A, A) == pytest.approx(0.0, abs=1e-5)


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
def test_multi_spyrmsd(minimize: bool, referenceRMSDs: List[float]) -> None:

    molc = copy.deepcopy(molecules.docking_1cbr[0])
    mols = [copy.deepcopy(mol) for mol in molecules.docking_1cbr[1:]]

    molc.strip()

    for mol in mols:
        mol.strip()

    RMSDs = rmsd.symmrmsd(
        molc.coordinates,
        [mol.coordinates for mol in mols],
        molc.atomicnums,
        mols[0].atomicnums,
        molc.adjacency_matrix,
        mols[0].adjacency_matrix,
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
def test_symmrmsd_cache(minimize: bool, referenceRMSDs: List[float]) -> None:

    molc = copy.deepcopy(molecules.docking_1cbr[0])
    mols = [copy.deepcopy(mol) for mol in molecules.docking_1cbr[1:]]

    molc.strip()

    for mol in mols:
        mol.strip()

    RMSDs = rmsd.symmrmsd(
        molc.coordinates,
        [mol.coordinates for mol in mols],
        molc.atomicnums,
        mols[0].atomicnums,
        molc.adjacency_matrix,
        mols[0].adjacency_matrix,
        minimize=minimize,
        cache=False,
    )

    for RMSD, referenceRMSD in zip(RMSDs, referenceRMSDs):
        assert RMSD == pytest.approx(referenceRMSD, abs=1e-5)


def test_issue_35_1():
    """
    GitHub Issue #35 from @kjelljorner

    https://github.com/RMeli/spyrmsd/issues/35
    """

    elements = np.array([6, 1, 1, 1, 1])

    connectivity_matrix = np.array(
        [
            [0, 1, 1, 1, 1],
            [1, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
        ]
    )

    coordinates_1 = np.array(
        [
            [2.416690358222e-09, -2.979634307864e-08, -9.782357562900e-09],
            [-6.669665490935e-01, -6.162099127861e-01, -6.069106091317e-01],
            [-4.258088724928e-01, 9.999808728518e-01, 1.078170393201e-01],
            [1.178459756093e-01, -4.538168967035e-01, 9.864390766805e-01],
            [9.749294435604e-01, 7.004596643409e-02, -4.873454970866e-01],
        ]
    )

    coordinates_2 = np.array(
        [
            [-2.118450971480e-07, 2.238951108509e-07, 1.839989120690e-07],
            [-5.297519571039e-01, -4.011375110922e-01, 8.668054003529e-01],
            [-5.107749001064e-01, 8.975573096842e-01, -3.555275589573e-01],
            [1.644944812511e-02, -7.486078704316e-01, -7.951194721576e-01],
            [1.024077620930e00, 2.521878479445e-01, 2.838414467631e-01],
        ]
    )

    r = rmsd.symmrmsd(
        coordinates_1,
        coordinates_2,
        elements,
        elements,
        connectivity_matrix,
        connectivity_matrix,
        center=True,
        minimize=True,
    )

    assert r == pytest.approx(0.0)


def test_issue_35_2():
    """
    GitHub Issue #35 from @kjelljorner

    https://github.com/RMeli/spyrmsd/issues/35
    """
    elements = np.array([6, 6, 1, 1, 1, 1, 1, 1])

    connectivity_matrix = np.array(
        [
            [0, 1, 1, 1, 1, 0, 0, 0],
            [1, 0, 0, 0, 0, 1, 1, 1],
            [1, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0],
        ]
    )

    coordinates_1 = np.array(
        [
            [-0.7513646148641476, 0.0212153090068301, 0.0811236469164399],
            [0.7513646329421021, -0.0212139873589844, -0.0811234544788992],
            [-1.1912097695762844, -0.9520810297276773, -0.1560519813078418],
            [-1.0197601491429782, 0.2790338755490671, 1.1099648032811764],
            [-1.1891251185534542, 0.769045487934088, -0.5868099187409093],
            [1.1891224079259752, -0.7690462686978631, 0.5868095572711365],
            [1.0197615786233758, -0.2790344505405097, -1.1099636171164908],
            [1.191211032645396, 0.9520810638350493, 0.1560509641753798],
        ]
    )

    coordinates_2 = np.array(
        [
            [0.752989230839992, 0.0675001892809206, -0.0055472904918074],
            [-0.7529892426611152, -0.0675008847728476, 0.0055489067798057],
            [1.0506759954349485, 1.0138756198991818, -0.4668239152871582],
            [1.2080371182091196, -0.7501873566674166, -0.5724157781772267],
            [1.148732484392816, 0.0417675099494674, 1.0141322717252037],
            [-1.050678518706957, -1.0138763219763327, 0.4668256362472835],
            [-1.148731454954795, -0.0417698443724041, -1.0141315673229787],
            [-1.208035612554004, 0.7501910886594293, 0.5724117365268774],
        ],
    )

    mask = elements != 1

    r = rmsd.symmrmsd(
        coordinates_1[mask],
        coordinates_2[mask],
        elements[mask],
        elements[mask],
        connectivity_matrix[mask, :][:, mask],
        connectivity_matrix[mask, :][:, mask],
        center=True,
        minimize=True,
    )

    assert r == pytest.approx(0.0)


@pytest.mark.parametrize(
    "i, j, result", [(1, 2, 1.95277757), (1, 3, 3.11801105), (2, 3, 2.98609758)]
)
def test_rmsd_atol(i: int, j: int, result: float):
    """
    Test usage of the :code:`atol` parameter for the QCP method.

    This parameter has been exposed to users following Issue 35 from @kjelljorner
    (https://github.com/RMeli/spyrmsd/issues/35)
    """

    moli = copy.deepcopy(molecules.docking_2viz[i])
    molj = copy.deepcopy(molecules.docking_2viz[j])

    # Check results are different from 0.0
    assert not result == pytest.approx(0.0)

    assert (
        rmsd.rmsd(
            moli.coordinates,
            molj.coordinates,
            moli.atomicnums,
            molj.atomicnums,
            minimize=True,
        )
        == pytest.approx(result)
    )

    assert (
        rmsd.rmsd(
            moli.coordinates,
            molj.coordinates,
            moli.atomicnums,
            molj.atomicnums,
            minimize=True,
            atol=1e9,
        )
        == pytest.approx(0.0)
    )


# Results obtained with OpenBabel
@pytest.mark.parametrize("i, reference", [(1, 0.476858), (2, 1.68089), (3, 1.50267)])
def test_symmrmsd_atol(i: bool, reference: float) -> None:

    moli = copy.deepcopy(molecules.docking_1cbr[0])
    molj = copy.deepcopy(molecules.docking_1cbr[i])

    moli.strip()
    molj.strip()

    # Check results are different from 0.0
    assert not reference == pytest.approx(0.0)

    assert (
        rmsd.symmrmsd(
            moli.coordinates,
            molj.coordinates,
            moli.atomicnums,
            molj.atomicnums,
            moli.adjacency_matrix,
            molj.adjacency_matrix,
            minimize=True,
        )
        == pytest.approx(reference, abs=1e-5)
    )

    assert (
        rmsd.symmrmsd(
            moli.coordinates,
            molj.coordinates,
            moli.atomicnums,
            molj.atomicnums,
            moli.adjacency_matrix,
            molj.adjacency_matrix,
            minimize=True,
            atol=1e9,
        )
        == pytest.approx(0.0)
    )


def test_symmrmsd_atol_multi() -> None:

    references = [0.476858, 1.68089, 1.50267]

    molc = copy.deepcopy(molecules.docking_1cbr[0])
    mols = [copy.deepcopy(mol) for mol in molecules.docking_1cbr[1:4]]

    molc.strip()

    for mol in mols:
        mol.strip()

    # Check results are different from 0.0
    assert not np.allclose(references, 0.0)

    RMSDs = rmsd.symmrmsd(
        molc.coordinates,
        [mol.coordinates for mol in mols],
        molc.atomicnums,
        mols[0].atomicnums,
        molc.adjacency_matrix,
        mols[0].adjacency_matrix,
        minimize=True,
    )

    for r, ref in zip(RMSDs, references):
        assert r == pytest.approx(ref, abs=1e-5)

    RMSDs = rmsd.symmrmsd(
        molc.coordinates,
        [mol.coordinates for mol in mols],
        molc.atomicnums,
        mols[0].atomicnums,
        molc.adjacency_matrix,
        mols[0].adjacency_matrix,
        minimize=True,
        atol=1e9,
    )

    for r in RMSDs:
        assert r == pytest.approx(0.0)


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
    "minimize, referenceRMSDs",
    [
        (
            False,  # No minimize: dummy RMSD
            [
                4.812480551076202,
                6.772045449820714,
                9.344911262612964,
                9.772939589989000,
                8.901837608843241,
            ],
        ),
        (
            True,  # Minimize: QCP
            [
                1.6578281551053196,
                1.7175638492348284,
                1.5946081072641485,
                2.1234944939308220,
                2.4894805175766606,
            ],
        ),
    ],
)
def test_rmsdwrapper_nosymm_protein(minimize: bool, referenceRMSDs: List[float]):

    mol0 = copy.deepcopy(molecules.trp[0])
    mols = [copy.deepcopy(mol) for mol in molecules.trp[1:]]

    RMSDs = rmsd.rmsdwrapper(mol0, mols, symmetry=False, minimize=minimize, strip=False)

    for RMSD, referenceRMSD in zip(RMSDs, referenceRMSDs):
        assert RMSD == pytest.approx(referenceRMSD)


@pytest.mark.parametrize(
    # Reference results obtained with OpenBabel
    "minimize, referenceRMSDs",
    [
        (
            True,  # Minimize: QCP + Isomorphism
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
        (
            False,  # No minimize: Isomorphism only
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
    ],
)
def test_rmsdwrapper_isomorphic(minimize: bool, referenceRMSDs: List[float]) -> None:

    molref = copy.deepcopy(molecules.docking_1cbr[0])
    mols = [copy.deepcopy(mol) for mol in molecules.docking_1cbr[1:]]

    RMSDs = rmsd.rmsdwrapper(molref, mols, minimize=minimize, strip=True)

    for RMSD, referenceRMSD in zip(RMSDs, referenceRMSDs):
        assert RMSD == pytest.approx(referenceRMSD, abs=1e-5)


@pytest.mark.parametrize(
    # Reference results obtained with OpenBabel
    "minimize, referenceRMSD",
    [(True, 0.476858), (False, 0.592256)],
)
def test_rmsdwrapper_single_molecule(minimize: bool, referenceRMSD: float) -> None:

    molref = copy.deepcopy(molecules.docking_1cbr[0])
    mols = copy.deepcopy(molecules.docking_1cbr[1])

    RMSD = rmsd.rmsdwrapper(molref, mols, minimize=minimize, strip=True)

    assert RMSD[0] == pytest.approx(referenceRMSD, abs=1e-5)
