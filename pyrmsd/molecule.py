from pyrmsd import utils, graph

import qcelemental as qcel
import numpy as np

try:
    # 3.0
    from openbabel import openbabel as ob
    from openbabel import pybel
except ImportError:
    # 2.0
    import openbabel as ob
    import pybel


class Molecule:
    def __init__(self, atomicnums, coordinates, adjacency_matrix=None):

        atomicnums = np.asarray(atomicnums, dtype=int)
        coordinates = np.asarray(coordinates, dtype=float)

        self.natoms = len(atomicnums)

        assert atomicnums.shape == (self.natoms,)
        assert coordinates.shape == (self.natoms, 3)

        self.atomicnums = atomicnums
        self.coordinates = coordinates

        self.stripped = np.all(atomicnums != 1)

        if adjacency_matrix is not None:
            self.adjacency_matrix = adjacency_matrix

        self.G = None

    def translate(self, vector):
        assert len(vector) == 3
        # TODO: More efficient (numpy vectorsation)
        for i, coord in enumerate(self.coordinates):
            self.coordinates[i] += vector

    def rotate(self, angle, axis, units="rad"):
        assert len(axis) == 3
        # TODO: More efficient (numpy vectorsation)
        for i, coord in enumerate(self.coordinates):
            self.coordinates[i] = utils.rotate(coord, angle, axis, units)

    def center_of_mass(self):
        # TODO: Cache this
        masses = [qcel.periodictable.to_mass(anum) for anum in self.atomicnums]

        return np.average(self.coordinates, axis=0, weights=masses)

    def center_of_geometry(self):
        return np.mean(self.coordinates, axis=0)

    def strip(self):

        if not self.stripped:
            idx = self.atomicnums != 1  # Hydrogen atoms

            # Strip
            self.atomicnums = self.atomicnums[idx]
            self.coordinates = self.coordinates[idx, :]

            # Update number of atoms
            self.natoms = len(self.atomicnums)

            self.stripped = True

    def to_graph(self):
        raise NotImplementedError

    def __len__(self):
        return self.natoms


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

    fmt = utils.format_openbabel(fname)

    obmol = next(pybel.readfile(fmt, fname))

    return obmol


def openbabel_to_molecule(obmol: ob.OBMol, adjacency: bool = True) -> Molecule:
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
        A = graph.adjacency_matrix_from_obmol(obmol)

    return Molecule(atomicnums, coordinates, A)
