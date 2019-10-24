import numpy as np
import scipy


def cost_mtx(A, B):

    return scipy.spatial.distance.cdist(A, B) ** 2


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
