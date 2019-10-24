import numpy as np

from scipy import optimize


def M_mtx(A, B):
    return B.T @ A


def K_mtx(M):

    assert M.shape == (3, 3)

    S_xx = M[0, 0]
    S_xy = M[0, 1]
    S_xz = M[0, 2]
    S_yx = M[1, 0]
    S_yy = M[1, 1]
    S_yz = M[1, 2]
    S_zx = M[2, 0]
    S_zy = M[2, 1]
    S_zz = M[2, 2]

    # p = plus, m = minus
    S_xx_yy_zz_ppp = S_xx + S_yy + S_zz
    S_yz_zy_pm = S_yz - S_zy
    S_zx_xz_pm = S_zx - S_xz
    S_xy_yx_pm = S_xy - S_yx
    S_xx_yy_zz_pmm = S_xx - S_yy - S_zz
    S_xy_yx_pp = S_xy + S_yx
    S_zx_xz_pp = S_zx + S_xz
    S_xx_yy_zz_mpm = -S_xx + S_yy - S_zz
    S_yz_zy_pp = S_yz + S_zy
    S_xx_yy_zz_mmp = -S_xx - S_yy + S_zz

    return np.array(
        [
            [S_xx_yy_zz_ppp, S_yz_zy_pm, S_zx_xz_pm, S_xy_yx_pm],
            [S_yz_zy_pm, S_xx_yy_zz_pmm, S_xy_yx_pp, S_zx_xz_pp],
            [S_zx_xz_pm, S_xy_yx_pp, S_xx_yy_zz_mpm, S_yz_zy_pp],
            [S_xy_yx_pm, S_zx_xz_pp, S_yz_zy_pp, S_xx_yy_zz_mmp],
        ]
    )


def coefficients(M, K):

    c2 = -2 * np.trace(M.T @ M)
    c1 = -8 * np.linalg.det(M)  # TODO: Slow?
    c0 = np.linalg.det(K)  # TODO: Slow?

    return c2, c1, c0


def lambda_max(Ga, Gb, c2, c1, c0):
    def P(x):
        return x ** 4 + c2 * x ** 2 + c1 * x + c0

    def dP(x):
        return 4 * x ** 3 + 2 * c2 * x + c1

    x0 = (Ga + Gb) / 2.0

    return optimize.newton(P, x0, fprime=dP)


def qcp_rmsd(A, B):

    assert A.shape == B.shape

    N = A.shape[0]

    Ga = np.trace(A.T @ A)
    Gb = np.trace(B.T @ B)

    M = M_mtx(A, B)
    K = K_mtx(M)

    assert np.allclose(K, K.T)

    c2, c1, c0 = coefficients(M, K)

    l_max = lambda_max(Ga, Gb, c2, c1, c0)

    s = Ga + Gb - 2 * l_max
    if abs(s) < 1e-12:  # Avoid numerical errors when Ga + Gb = 2 * l_max
        rmsd = 0
    else:
        rmsd = np.sqrt(s / N)

    return rmsd
