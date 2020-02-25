from typing import List, Tuple

import numpy as np
import rdkit.Chem as Chem

from spyrmsd import molecule, utils


def load(fname: str):
    """
    Load molecule from file.

    Parameters
    ----------
    fname: str
        File name

    Returns
    -------
    Molecule
    """

    fmt = utils.molformat(fname)

    if fmt == "mol2":
        rdmol = Chem.MolFromMol2File(fname, removeHs=False)
    elif fmt == "sdf":
        rdmol = Chem.SDMolSupplier(fname, removeHs=False)[0]
    elif fmt == "pdb":
        rdmol = Chem.MolFromPDBFile(fname, removeHs=False)
    else:
        raise NotImplementedError

    return rdmol


def loadall(fname: str):
    """
    Load molecules from file.

    Parameters
    ----------
    fname: str
        File name

    Returns
    -------
    List of molecules
    """

    fmt = utils.molformat(fname)

    if fmt == "mol2":
        raise NotImplementedError  # See RDKit Issue #415
    elif fmt == "sdf":
        rdmols = Chem.SDMolSupplier(fname, removeHs=False)
    elif fmt == "pdb":
        # TODO: Implement
        raise NotImplementedError
    else:
        raise NotImplementedError

    return [rdmol for rdmol in rdmols]


def adjacency_matrix(mol) -> np.ndarray:
    """
    Adjacency matrix from OpenBabel molecule.

    Parameters
    ----------
    mol:
        Molecule

    Returns
    -------
    np.ndarray
        Adjacency matrix of the molecule
    """

    return Chem.rdmolops.GetAdjacencyMatrix(mol)


def to_molecule(mol, adjacency: bool = True) -> molecule.Molecule:
    """
    Transform molecule to `pyrmsd` molecule.

    Parameters
    ----------
    mol:
        Molecule
    adjacency: boolean, optional
        Flag to decide wether to build the adjacency matrix from molecule

    Returns
    -------
    spyrmsd.molecule.Molecule
        `spyrmsd` molecule
    """

    atoms = mol.GetAtoms()

    n = len(atoms)

    atomicnums = np.zeros(n, dtype=int)
    coordinates = np.zeros((n, 3))

    conformer = mol.GetConformer()

    for i, atom in enumerate(atoms):
        atomicnums[i] = atom.GetAtomicNum()

        pos = conformer.GetAtomPosition(i)

        coordinates[i] = np.array([pos.x, pos.y, pos.z])

    if adjacency:
        A = adjacency_matrix(mol)

    return molecule.Molecule(atomicnums, coordinates, A)


def numatoms(mol) -> int:
    """
    Number of atoms.

    Parameters
    ----------
    mol:
        Molecule

    Returns
    -------
    int
        Number of atoms
    """
    return mol.GetNumAtoms()


def numbonds(mol) -> int:
    """
    Number of bonds.

    Parameters
    ----------
    mol:
        Molecule

    Returns
    -------
    int
        Number of bonds
    """
    return mol.GetNumBonds()


def bonds(mol) -> List[Tuple[int, int]]:
    """
    List of bonds.

    Parameters
    ----------
    mol:
        Molecule

    Returns
    -------
    List[Tuple[int, int]]
        List of bonds

    Notes
    -----
    A bond is defined by a tuple of (0-based) indices of two atoms.
    """
    b = []

    for bond in mol.GetBonds():
        b.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

    return b
