from pyrmsd import molecule, io

import os

try:
    from openbabel import pybel  # 3.0
except ImportError:
    import pybel  # 2.0

from typing import Tuple

fdir = os.path.dirname(os.path.abspath(__file__))
molpath = os.path.join(fdir, "data/molecules/")


def load(fname: str) -> Tuple[pybel.Molecule, molecule.Molecule]:
    """
    Load molecule from file.

    Parameters
    ----------
    fname: str
        Input file name

    Returns
    -------
    Tuple[pybel.Molecule, molecule.Molecule]
        Loaded molecule as `pybel.Molecule` and `pyrmsd.molecule.Molecule`
    """

    fname = os.path.join(molpath, fname)

    obmol = io.load(fname)

    mol = io.openbabel_to_molecule(obmol, adjacency=True)

    return obmol, mol


obbenzene, benzene = load("benzene.xyz")
obethanol, ethanol = load("ethanol.xyz")
obxyz = [obbenzene, obethanol]
xyz = [benzene, ethanol]

obdialanine, dialanine = load("dialanine.sdf")
obsdf = [obdialanine]
sdf = [dialanine]

allmolecules = xyz + sdf
allobmolecules = obxyz + obsdf

docking_2viz = {}
obdocking_2viz = {}
for i in [1, 2, 3]:
    obdocking_2viz[i], docking_2viz[i] = load(f"2viz_{i}.sdf")
