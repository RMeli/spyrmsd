import os
from typing import Any, List, Tuple

from spyrmsd import io, molecule

fdir = os.path.dirname(os.path.abspath(__file__))
molpath = os.path.join(fdir, "data/molecules/")


def load(fname: str) -> Tuple[Any, molecule.Molecule]:
    """
    Load molecule from file.

    Parameters
    ----------
    fname: str
        Input file name

    Returns
    -------
    Tuple[Any, molecule.Molecule]
        Loaded molecule as `pybel.Molecule` or `rdkit.Chem.rdkem.Mol` and
        `pyrmsd.molecule.Molecule`
    """

    fname = os.path.join(molpath, fname)

    m = io.load(fname)

    mol = io.to_molecule(m, adjacency=True)

    return m, mol


def loadall(fname: str) -> Tuple[List[Any], List[molecule.Molecule]]:
    """
    Load all molecule from file.

    Parameters
    ----------
    fname: str
        Input file name

    Returns
    -------
    Tuple[List[Any], List[molecule.Molecule]]
        Loaded molecule as `pybel.Molecule` or `rdkit.Chem.rdchem.Mol` and
        `pyrmsd.molecule.Molecule`
    """

    fname = os.path.join(molpath, fname)

    ms = io.loadall(fname)

    mols = [io.to_molecule(m, adjacency=True) for m in ms]

    return ms, mols


obbenzene, benzene = load("benzene.sdf")
obpyridine, pyridine = load("pyridine.sdf")
obethanol, ethanol = load("ethanol.sdf")
obdialanine, dialanine = load("dialanine.sdf")
obsdf = [obbenzene, obpyridine, obethanol, obdialanine]
sdf = [benzene, pyridine, ethanol, dialanine]

allmolecules = sdf
allobmolecules = obsdf

obdocking_2viz, docking_2viz = {}, {}
for i in [1, 2, 3]:
    obdocking_2viz[i], docking_2viz[i] = load(f"2viz_{i}.sdf")

obdocking_1cbr = [load("1cbr_ligand.mol2")[0], *loadall("1cbr_docking.sdf")[0]]
docking_1cbr = [load("1cbr_ligand.mol2")[1], *loadall("1cbr_docking.sdf")[1]]

intrp, trp = [], []
for i in range(6):
    intrp.append(load(f"trp{i}.pdb")[0])
    trp.append(load(f"trp{i}.pdb")[1])
