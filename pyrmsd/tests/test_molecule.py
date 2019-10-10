from pyrmsd import molecule

import os

fdir = os.path.dirname(os.path.abspath(__file__))
molpath = os.path.join(fdir, os.pardir, "data/molecules/")

def load(fname: str):

    fname = os.path.join(molpath, fname)

    mol = molecule.load(fname)

    return mol

benzene = load("benzene.xyz")

def test_load_benzene():

    mol = benzene

    assert len(mol.atoms) == 12

    num_h = num_c = 0
    for atom in mol.atoms:
        if atom.atomicnum == 1:
            num_h += 1
        if atom.atomicnum == 6:
            num_c += 1
    assert num_h == num_c == 6

def test_openbabel_to_molecule_benzene():

    mol = benzene

    m = molecule.openbabel_to_molecule(mol)

    assert m.atomicnums.shape == (12,)
    assert m.coordinates.shape == (12, 3)

