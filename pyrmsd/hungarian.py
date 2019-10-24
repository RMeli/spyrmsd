import numpy as np
import scipy


def cost_mtx(A, B):

    assert A.shape == B.shape
    assert A.shape[1] == B.shape[1] == 3

    # Vectorization of pairwise distances
    # (B-A)**2 = A**2  + B**2 - 2*A*B
    A2 = np.sum(A ** 2, axis=1)
    B2 = np.sum(B ** 2, axis=1)[:, np.newaxis]
    C = -2 * A @ B.T + A2 + B2

    # Subtraction can lead to small negative elements because of numerical errors
    C[np.abs(C) < 1e-12] = 0.0

    return C


def optimal_assignment(A, B):

    C = cost_mtx(A, B)

    row_idx, col_idx = scipy.optimize.linear_sum_assignment(C)

    # Compute assignment cost
    cost = C[row_idx, col_idx].sum()

    return cost, row_idx, col_idx


def hungarian_rmsd(A, B, typesA, typesB):

    assert A.shape == B.shape
    assert typesA.shape == typesB.shape

    types = set(typesA)

    total_cost: float = 0.0
    for t in types:
        typesA_idx = typesA == t
        typesB_idx = typesB == t

        assert typesA_idx.shape == typesB_idx.shape

        cost, row_idx, col_idx = optimal_assignment(A[typesA_idx, :], B[typesB_idx, :])

        total_cost += cost

    N = A.shape[0]

    rmsd = np.sqrt(total_cost / N)

    return rmsd
