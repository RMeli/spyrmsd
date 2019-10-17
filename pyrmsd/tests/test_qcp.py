from pyrmsd import rmsd, molecule, qcp
from pyrmsd.tests import molecules

import numpy as np

import pytest


def test_qcm_benzene():

    mol1 = molecules.benzene
    mol2 = molecules.benzene

    m1 = molecule.openbabel_to_molecule(mol1)
    m2 = molecule.openbabel_to_molecule(mol2)

    assert rmsd.rmsd_dummy(m1, m2) == pytest.approx(0)

    m2.rotate(10, np.array([1, -2, -1]))

    assert rmsd.rmsd_dummy(m1, m2) > 0.0

    assert qcp.qcp_rmsd(m1.coordinates, m2.coordinates) == pytest.approx(0)
