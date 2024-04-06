import copy
import itertools
from typing import Tuple

import numpy as np
import pytest

from spyrmsd import qcp


def test_M_mtx(mol) -> None:
    mol1 = copy.deepcopy(mol.mol)
    mol2 = copy.deepcopy(mol.mol)

    # Build rotated coordinate set
    mol2.rotate(10, np.random.rand(3))

    M = qcp.M_mtx(mol1.coordinates, mol2.coordinates)

    assert M.shape == (3, 3)

    def S(i, j):
        A, B = mol1.coordinates, mol2.coordinates
        return np.dot(B[:, i], A[:, j])

    for i, j in itertools.combinations(range(3), 2):
        assert M[i, j] == pytest.approx(S(i, j))


def test_K_mtx(mol) -> None:
    mol1 = copy.deepcopy(mol.mol)
    mol2 = copy.deepcopy(mol.mol)

    # Build rotated coordinate set
    mol2.rotate(10, np.random.rand(3))

    M = qcp.M_mtx(mol1.coordinates, mol2.coordinates)
    K = qcp.K_mtx(M)

    assert K.shape == (4, 4)
    assert np.allclose(K.T, K)
    assert np.trace(K) == pytest.approx(0)


@pytest.mark.parametrize(
    "input, result",
    [
        ((2, 2, 0, 0, -1), 1),  # f(x) = x^4 - 1; x_0 = 2
        ((-2, -2, 0, 0, -1), -1),  # f(x) = x^4 - 1; x_0 = -2
        ((3, 3, -5, 0, 4), 2),  # f(x) = x^4 - 5 * x^2 + 4; x_0 = 3
        ((1, 1, -5, 0, 4), 1),  # f(x) = x^4 - 5 * x^2 + 4; x_0 = 1/2
        ((-1, -1, -5, 0, 4), -1),  # f(x) = x^4 - 5 * x^2 + 4; x_0 = -1/2
        ((-3, -3, -5, 0, 4), -2),  # f(x) = x^4 - 5 * x^2 + 4; x_0 = -3
    ],
    ids=["f1", "f2", "f3", "f4", "f5", "f6"],
)
def test_lambda_max(
    input: Tuple[float, float, float, float, float], result: float
) -> None:
    assert qcp.lambda_max(*input) == pytest.approx(result)


def test_lambda_max_eig(mol) -> None:
    mol1 = copy.deepcopy(mol.mol)
    mol2 = copy.deepcopy(mol.mol)

    # Build rotated coordinate set
    mol2.rotate(10, np.random.rand(3))

    A = mol1.coordinates
    B = mol2.coordinates

    Ga = np.trace(A.T @ A)
    Gb = np.trace(B.T @ B)

    M = qcp.M_mtx(A, B)
    K = qcp.K_mtx(M)

    assert K.shape == (4, 4)
    assert np.allclose(K.T, K)
    assert np.trace(K) == pytest.approx(0)

    c2, c1, c0 = qcp.coefficients(M, K)

    lm = qcp.lambda_max(Ga, Gb, c2, c1, c0)
    lm_eig = qcp._lambda_max_eig(K)

    assert lm_eig == pytest.approx(lm)
