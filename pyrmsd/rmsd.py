import numpy as np


def rmsd_dummy(mol1, mol2, center=False):

    assert np.all(mol1.atomicnums == mol2.atomicnums)

    n = len(mol1)

    c1 = mol1.coordinates
    c2 = mol2.coordinates

    if center:
        c1 -= mol1.center_of_geometry()
        c2 -= mol2.center_of_geometry()

    return np.sqrt(np.sum((c1 - c2) ** 2) / n)
