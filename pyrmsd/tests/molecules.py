from pyrmsd import molecule

import os

fdir = os.path.dirname(os.path.abspath(__file__))
molpath = os.path.join(fdir, os.pardir, "data/molecules/")

def load(fname: str):

    fname = os.path.join(molpath, fname)

    mol = molecule.load(fname)

    return mol

benzene = load("benzene.xyz")