from pyrmsd import molecule

import os

fdir = os.path.dirname(os.path.abspath(__file__))
molpath = os.path.join(fdir, os.pardir, "data/molecules/")


def load(fname: str) -> molecule.Molecule:

    fname = os.path.join(molpath, fname)

    m = molecule.load(fname)

    mol = molecule.openbabel_to_molecule(m)

    return mol


benzene = load("benzene.xyz")
ethanol = load("ethanol.xyz")
xyz = [benzene, ethanol]

dialanine = load("dialanine.sdf")
sdf = [dialanine]

docking_2viz = {}
for i in [1, 2, 3]:
    docking_2viz[i] = load(f"2viz_{i}.sdf")

allmolecules = xyz + sdf
