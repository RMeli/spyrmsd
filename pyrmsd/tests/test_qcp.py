from pyrmsd import molecule, qcp
from pyrmsd.tests import molecules

import numpy as np
import itertools

import pytest


def test_M_mtx():
    for mol in molecules.allmolecules:
        mol1 = molecule.openbabel_to_molecule(mol)
        mol2 = molecule.openbabel_to_molecule(mol)

        # Build rotated coordinate set
        mol2.rotate(10, np.random.rand(3))

        M = qcp.M_mtx(mol1.coordinates, mol2.coordinates)

        assert M.shape == (3, 3)

        def S(i, j):
            A, B = mol1.coordinates, mol2.coordinates
            return np.dot(B[:, i], A[:, j])

        for i, j in itertools.combinations(range(3), 2):
            assert M[i, j] == pytest.approx(S(i, j))


def test_K_mtx():

    for mol in molecules.allmolecules:
        mol1 = molecule.openbabel_to_molecule(mol)
        mol2 = molecule.openbabel_to_molecule(mol)

        # Build rotated coordinate set
        mol2.rotate(10, np.random.rand(3))

        M = qcp.M_mtx(mol1.coordinates, mol2.coordinates)
        K = qcp.K_mtx(M)

        assert K.shape == (4, 4)
        assert np.allclose(K.T, K)
        assert np.trace(K) == pytest.approx(0)


def test_lambda_max():

    # f(x) = x^4 - 1; x_0 = 2
    assert qcp.lambda_max(2, 2, 0, 0, -1) == pytest.approx(1)

    # f(x) = x^4 - 1; x_0 = -2
    assert qcp.lambda_max(-2, -2, 0, 0, -1) == pytest.approx(-1)

    # f(x) = x^4 - 5 * x^2 + 4; x_0 = 3
    assert qcp.lambda_max(3, 3, -5, 0, 4) == pytest.approx(2)

    # f(x) = x^4 - 5 * x^2 + 4; x_0 = 1/2
    assert qcp.lambda_max(1, 1, -5, 0, 4) == pytest.approx(1)

    # f(x) = x^4 - 5 * x^2 + 4; x_0 = -1/2
    assert qcp.lambda_max(-1, -1, -5, 0, 4) == pytest.approx(-1)

    # f(x) = x^4 - 5 * x^2 + 4; x_0 = -3
    assert qcp.lambda_max(-3, -3, -5, 0, 4) == pytest.approx(-2)
