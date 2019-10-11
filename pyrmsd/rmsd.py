import numpy as np


def rmsd_molecules(mol1, mol2):

    n = len(mol1)

    assert np.all(mol1.atomicnums == mol2.atomicnums)

    c1 = mol1.coordinates
    c2 = mol2.coordinates

    return np.sqrt(np.sum((c1 - c2) ** 2) / n)
