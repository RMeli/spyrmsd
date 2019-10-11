from pyrmsd import rmsd, molecule
from pyrmsd.tests import molecules

import pytest


def test_rmsd_benzene():

    mol1 = molecules.benzene
    mol2 = molecules.benzene

    m1 = molecule.openbabel_to_molecule(mol1)
    m2 = molecule.openbabel_to_molecule(mol2)

    assert rmsd.rmsd_molecules(m1, m2) == pytest.approx(0)
