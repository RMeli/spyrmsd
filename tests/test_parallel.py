import copy
from typing import List

import numpy as np
import pytest

try:
    from spyrmsd.parallel import prmsdwrapper
except ImportError:
    pytest.skip("Pebble not available", allow_module_level=True)


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
    ids=["no_minimize", "minimize"],
)
def test_prmsdwrapper_nosymm_protein(trps, minimize: bool, referenceRMSDs: List[float]):
    mol0 = copy.deepcopy(trps[0])
    mols = [copy.deepcopy(mol) for mol in trps[1:]]

    RMSDs = prmsdwrapper(mol0, mols, symmetry=False, minimize=minimize, strip=False)

    for RMSD, referenceRMSD in zip(RMSDs, referenceRMSDs):
        np.testing.assert_allclose(RMSD, referenceRMSD, atol=1e-5)


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
    ids=["minimize", "no_minimize"],
)
def test_prmsdwrapper_isomorphic(
    docking_1cbr, minimize: bool, referenceRMSDs: List[float]
) -> None:
    molref = copy.deepcopy(docking_1cbr[0])
    mols = [copy.deepcopy(mol) for mol in docking_1cbr[1:]]

    RMSDs = prmsdwrapper(molref, mols, minimize=minimize, strip=True)

    for RMSD, referenceRMSD in zip(RMSDs, referenceRMSDs):
        np.testing.assert_allclose(RMSD, referenceRMSD, atol=1e-5)


@pytest.mark.parametrize(
    # Reference results obtained with OpenBabel
    "minimize, referenceRMSD",
    [(True, 0.476858), (False, 0.592256)],
    ids=["minimize", "no_minimize"],
)
def test_prmsdwrapper_single_molecule(
    docking_1cbr, minimize: bool, referenceRMSD: float
) -> None:
    molref = copy.deepcopy(docking_1cbr[0])
    mols = copy.deepcopy(docking_1cbr[1])

    RMSD = prmsdwrapper(molref, mols, minimize=minimize, strip=True)

    np.testing.assert_allclose(RMSD[0], referenceRMSD, atol=1e-5)


def test_prmsdwrapper_single_molecule_timeout(muparfostat) -> None:
    mol1 = muparfostat.mol
    mol2 = muparfostat.mol

    with pytest.warns(
        UserWarning, match="1 chunks timed out"
    ):
        RMSD = prmsdwrapper(mol1, mol2, timeout=1e-3, num_workers=1)

    assert np.isnan(RMSD[0])


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
    ids=["minimize", "no_minimize"],
)
def test_prmsdwrapper_molecules_chunksize_no_timeout(
    docking_1cbr, minimize: bool, referenceRMSDs: List[float]
) -> None:
    molref = copy.deepcopy(docking_1cbr[0])
    mols = [copy.deepcopy(mol) for mol in docking_1cbr[1:]]

    RMSDlist = prmsdwrapper(molref, mols, minimize=minimize, chunksize=4, num_workers=1)

    for RMSD, referenceRMSD in zip(RMSDlist, referenceRMSDs):
        np.testing.assert_allclose(RMSD, referenceRMSD, atol=1e-5)


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
    ids=["minimize", "no_minimize"],
)
def test_prmsdwrapper_molecules_chunksize_timeout(
    docking_1cbr, minimize: bool, referenceRMSDs: List[float]
) -> None:
    molref = copy.deepcopy(docking_1cbr[0])
    mols = [copy.deepcopy(mol) for mol in docking_1cbr[1:]]

    with pytest.warns(UserWarning, match="When using the timeout feature"):
        RMSDlist = prmsdwrapper(
            molref,
            mols,
            minimize=minimize,
            timeout=1,
            chunksize=4,
            num_workers=1,
        )

    for RMSD, referenceRMSD in zip(RMSDlist, referenceRMSDs):
        np.testing.assert_allclose(RMSD, referenceRMSD, atol=1e-5)
