from pyrmsd import molecule, utils

from typing import List

import numpy as np

try:
    # 3.0
    from openbabel import openbabel as ob
    from openbabel import pybel
except ImportError:
    # 2.0
    import openbabel as ob
    import pybel


def load(fname: str) -> ob.OBMol:
    """
    Load molecule from file

    Parameters
    ----------
    fname: str
        File name

    Returns
    -------
    openbabel.OBMol
        OpenBabel molecule
    """

    fmt = utils.molformat(fname)

    obmol = next(pybel.readfile(fmt, fname))

    return obmol


def loadall(fname: str) -> List[ob.OBMol]:
    """
    Load molecules from file

    Parameters
    ----------
    fname: str
        File name

    Returns
    -------
    List[openbabel.OBMol]
        List of OpenBabel molecules
    """

    fmt = utils.molformat(fname)

    return [obmol for obmol in pybel.readfile(fmt, fname)]


def adjacency_matrix(obmol: ob.OBMol) -> np.ndarray:
    """
    Adjacency matrix from OpenBabel molecule.

    Parameters
    ----------
    openbabel.OBMol
        OpenBabel molecule

    Returns
    -------
    np.ndarray
        Adjacency matrix of the molecule
    """

    n = len(obmol.atoms)

    # Pre-allocate memory for  the adjacency matrix
    A = np.zeros((n, n), dtype=int)

    # Loop over molecular bonds
    for bond in ob.OBMolBondIter(obmol.OBMol):
        # Bonds are 1-indexed
        i: int = bond.GetBeginAtomIdx() - 1
        j: int = bond.GetEndAtomIdx() - 1

        # A molecular graph is undirected
        A[i, j] = A[j, i] = 1

    return A


def to_molecule(obmol: ob.OBMol, adjacency: bool = True) -> molecule.Molecule:
    """
    Transform OpenBabel molecule to `pyrmsd` molecule

    Parameters
    ----------
    obmol: ob.OBMol
        OpenBabel molecule
    adjacency: boolean, optional
        Flag to decide wether to build the adjacency matrix from the OpenBabel molecule

    Returns
    -------
    pyrmsd.molecule.Molecule
        `pyrmsd` molecule
    """

    n = len(obmol.atoms)

    atomicnums = np.zeros(n, dtype=int)
    coordinates = np.zeros((n, 3))

    for i, atom in enumerate(obmol.atoms):
        atomicnums[i] = atom.atomicnum
        coordinates[i] = atom.coords

    if adjacency:
        A = adjacency_matrix(obmol)

    return molecule.Molecule(atomicnums, coordinates, A)
