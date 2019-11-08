from pyrmsd import molecule, io

import os

try:
    from openbabel import pybel  # 3.0
except ImportError:
    import pybel  # 2.0

from typing import Tuple, List

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


def loadall(fname: str) -> Tuple[List[pybel.Molecule], List[molecule.Molecule]]:
    """
    Load all molecule from file.

    Parameters
    ----------
    fname: str
        Input file name

    Returns
    -------
    Tuple[List[pybel.Molecule], List[molecule.Molecule]]
        Loaded molecule as `pybel.Molecule` and `pyrmsd.molecule.Molecule`
    """

    fname = os.path.join(molpath, fname)

    obmols = io.loadall(fname)

    mols = [io.openbabel_to_molecule(obmol, adjacency=True) for obmol in obmols]

    return obmols, mols


obbenzene, benzene = load("benzene.xyz")
obethanol, ethanol = load("ethanol.xyz")
obxyz = [obbenzene, obethanol]
xyz = [benzene, ethanol]

obdialanine, dialanine = load("dialanine.sdf")
obsdf = [obdialanine]
sdf = [dialanine]

allmolecules = xyz + sdf
allobmolecules = obxyz + obsdf

obdocking_2viz, docking_2viz = {}, {}
for i in [1, 2, 3]:
    obdocking_2viz[i], docking_2viz[i] = load(f"2viz_{i}.sdf")

obdocking_1cbr = [load("1cbr_ligand.mol2")[0], *loadall("1cbr_docking.sdf")[0]]
docking_1cbr = [load("1cbr_ligand.mol2")[1], *loadall("1cbr_docking.sdf")[1]]

obtrp, trp = {}, {}
for i in range(6):
    obtrp[i], trp[i] = load(f"trp{i}.pdb")
