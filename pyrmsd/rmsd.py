import numpy as np


def rmsd_molecules(mol1, mol2):

    n = len(mol1)

    assert np.all(mol1.atomicnums == mol2.atomicnums)

    return np.sqrt(np.sum((mol1.coordinates - mol2.coordinates) ** 2) / n)
