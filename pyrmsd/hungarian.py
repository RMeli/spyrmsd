import numpy as np
import scipy


def cost_mtx(A, B):

    assert A.shape == B.shape
    assert A.shape[1] == B.shape[1] == 3

    n = A.shape[0]

    C = np.zeros((n, n))

    # TODO: Vectorize
    for i in range(n):
        for j in range(n):
            a = A[i, :]
            b = B[j, :]

            ab = b - a

            C[i, j] = np.dot(ab, ab)

    return C


def optimal_assignment(A, B):

    C = cost_mtx(A, B)

    row_idx, col_idx = scipy.optimize.linear_sum_assignment(C)

    # Compute assignment cost
    cost = C[row_idx, col_idx].sum()

    return cost, row_idx, col_idx


def hungarian_rmsd(A, B):

    assert A.shape == B.shape

    N = A.shape[0]

    cost, row_idx, col_idx = optimal_assignment(A, B)

    rmsd = np.sqrt(cost / N)

    return rmsd, row_idx, col_idx
