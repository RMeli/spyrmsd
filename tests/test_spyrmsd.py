import copy
import sys
from typing import List

import pytest

from spyrmsd import spyrmsd
from tests import molecules


def test_spyrmsd_imported():
    assert "spyrmsd" in sys.modules


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

    RMSDs = spyrmsd.rmsdwrapper(
        mol0, mols, symmetry=False, minimize=minimize, strip=False
    )

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

    RMSDs = spyrmsd.rmsdwrapper(molref, mols, minimize=minimize, strip=True)

    for RMSD, referenceRMSD in zip(RMSDs, referenceRMSDs):
        assert RMSD == pytest.approx(referenceRMSD, abs=1e-5)
