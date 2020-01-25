import copy

import numpy as np
import pytest

from spyrmsd import hungarian, molecule
from tests import molecules


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_cost_mtx(mol: molecule.Molecule):

    mol1 = copy.deepcopy(mol)
    mol2 = copy.deepcopy(mol)

    M = hungarian.cost_mtx(mol1.coordinates, mol2.coordinates)

    n = len(mol1)

    assert len(mol2) == n
    assert M.shape == (n, n)
    assert np.all(M >= 0)

    for i in range(n):
        for j in range(n):
            a, b = mol1.coordinates[i, :], mol2.coordinates[j, :]
            ab = b - a
            assert M[i, j] == pytest.approx(np.dot(ab, ab))


@pytest.mark.parametrize("mol", molecules.allmolecules)
def test_optimal_assignement_same_molecule(mol: molecule.Molecule):

    mol1 = copy.deepcopy(mol)
    mol2 = copy.deepcopy(mol)

    assert len(mol1) == len(mol2)

    n = len(mol1)

    cost, row_idx, col_idx = hungarian.optimal_assignment(
        mol1.coordinates, mol2.coordinates
    )

    # When molecules are the same, the assignment is perfect and the cost is zero
    assert cost == pytest.approx(0)

    assert len(row_idx) == len(col_idx) == n

    # The perfect assignment corresponds to the original ordering
    assert np.all(row_idx == np.arange(0, n))
    assert np.all(col_idx == np.arange(0, n))
