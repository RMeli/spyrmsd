from pyrmsd import qcp, hungarian, graph

import numpy as np


def rmsd_dummy(mol1, mol2, center=False):

    assert np.all(mol1.atomicnums == mol2.atomicnums)
    assert mol1.coordinates.shape == mol2.coordinates.shape

    n = len(mol1)

    c1 = mol1.coordinates
    c2 = mol2.coordinates

    if center:
        c1 -= mol1.center_of_geometry()
        c2 -= mol2.center_of_geometry()

    return np.sqrt(np.sum((c1 - c2) ** 2) / n)


def rmsd_qcp(mol1, mol2):

    assert np.all(mol1.atomicnums == mol2.atomicnums)

    c1 = mol1.coordinates - mol1.center_of_geometry()
    c2 = mol2.coordinates - mol2.center_of_geometry()

    return qcp.qcp_rmsd(c1, c2)


def rmsd_hungarian(mol1, mol2, center=False):

    assert mol1.atomicnums.shape == mol2.atomicnums.shape
    assert mol1.coordinates.shape == mol2.coordinates.shape

    c1 = mol1.coordinates
    c2 = mol2.coordinates

    if center:
        c1 -= mol1.center_of_geometry()
        c2 -= mol2.center_of_geometry()

    return hungarian.hungarian_rmsd(c1, c2, mol1.atomicnums, mol2.atomicnums)


def rmsd_isomorphic(mol1, mol2, center=False):

    assert mol1.atomicnums.shape == mol2.atomicnums.shape
    assert mol1.coordinates.shape == mol2.coordinates.shape

    n = len(mol1)

    c1 = mol1.coordinates
    c2 = mol2.coordinates

    if center:
        c1 -= mol1.center_of_geometry()
        c2 -= mol2.center_of_geometry()

    G1 = mol1.to_graph()
    G2 = mol2.to_graph()

    isomorphisms = graph.match_graphs(G1, G2)

    # Minimum squared displacement
    min_sd = np.inf
    for isomorphism in isomorphisms:

        # Use the isomorphism to shuffle coordinates around (from original order)
        c1i = c1[list(isomorphism.keys()), :]
        c2i = c2[list(isomorphism.values()), :]

        # Compute square displacement
        # Avoid dividing by n and an expensive sqrt() operation
        sd = np.sum((c1i - c2i) ** 2)

        if sd < min_sd:
            min_sd = sd

    return np.sqrt(min_sd / n)
