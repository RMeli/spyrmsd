import copy
import sys

import pytest

from spyrmsd import molecule, spyrmsd
from tests import molecules


def test_spyrmsd_imported():
    assert "spyrmsd" in sys.modules


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_rmsdwrapper_molsize(mol: molecule.Molecule) -> None:

    m = copy.deepcopy(mol)
    ms = copy.deepcopy(mol)

    ms.strip()

    with pytest.raises(ValueError):
        spyrmsd.rmsdwrapper(m, ms, symmetry=False)


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
def test_rmsdwrapper_qcp_protein(i: int, rmsd_dummy: float, rmsd_min: float):

    mol0 = copy.deepcopy(molecules.trp[0])
    mol = copy.deepcopy(molecules.trp[i])

    assert spyrmsd.rmsdwrapper(mol0, mol, symmetry=False) == pytest.approx(rmsd_dummy)
    assert spyrmsd.rmsdwrapper(
        mol0, mol, symmetry=False, minimize=True
    ) == pytest.approx(rmsd_min)


# Results obtained with OpenBabel
@pytest.mark.parametrize(
    "index, RMSD",
    [
        (1, 0.592256),
        (2, 2.11545),
        (3, 2.29824),
        (4, 9.45773),
        (5, 1.35005),
        (6, 9.44356),
        (7, 9.59758),
        (8, 9.55076),
        (9, 2.44067),
        (10, 9.6171),
    ],
)
def test_rmsdwrapper_isomorphic(index: int, RMSD: float) -> None:

    molc = copy.deepcopy(molecules.docking_1cbr[0])
    mol = copy.deepcopy(molecules.docking_1cbr[index])

    assert spyrmsd.rmsdwrapper(mol, molc, strip=True) == pytest.approx(RMSD, abs=1e-5)


# Results obtained with OpenBabel
@pytest.mark.parametrize(
    "index, RMSD",
    [
        (1, 0.476858),
        (2, 1.68089),
        (3, 1.50267),
        (4, 1.90623),
        (5, 1.01324),
        (6, 1.31716),
        (7, 1.11312),
        (8, 1.06044),
        (9, 0.965387),
        (10, 1.37842),
    ],
)
def test_rmsdwrapper_qcp_isomorphic(index: int, RMSD: float) -> None:

    molc = copy.deepcopy(molecules.docking_1cbr[0])
    mol = copy.deepcopy(molecules.docking_1cbr[index])

    assert spyrmsd.rmsdwrapper(mol, molc, minimize=True, strip=True) == pytest.approx(
        RMSD, abs=1e-5
    )
