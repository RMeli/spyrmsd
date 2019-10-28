from pyrmsd import molecule

import os

try:
    from openbabel import pybel  # 3.0
except ImportError:
    import pybel  # 2.0

from typing import Tuple

fdir = os.path.dirname(os.path.abspath(__file__))
molpath = os.path.join(fdir, os.pardir, "data/molecules/")


def load(fname: str) -> Tuple[pybel.Molecule, molecule.Molecule]:

    fname = os.path.join(molpath, fname)

    obmol = molecule.load(fname)

    mol = molecule.openbabel_to_molecule(obmol)

    return obmol, mol


obbenzene, benzene = load("benzene.xyz")
obethanol, ethanol = load("ethanol.xyz")
obxyz = [obbenzene, obethanol]
xyz = [benzene, ethanol]

obdialanine, dialanine = load("dialanine.sdf")
obsdf = [obdialanine]
sdf = [dialanine]

docking_2viz = {}
obdocking_2viz = {}
for i in [1, 2, 3]:
    obdocking_2viz[i], docking_2viz[i] = load(f"2viz_{i}.sdf")

allobmolecules = obxyz + obsdf + list(obdocking_2viz.values())
allmolecules = xyz + sdf
