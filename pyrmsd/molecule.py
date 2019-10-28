from pyrmsd import utils

import qcelemental as qcel
import numpy as np

try:
    import pybel  # 2.0
except ImportError:
    from openbabel import pybel  # 3.0


def load(fname: str):

    fmt = utils.format_openbabel(fname)

    obmol = next(pybel.readfile(fmt, fname))

    return obmol


def openbabel_to_molecule(obmol):

    n = len(obmol.atoms)

    atomicnums = np.zeros(n, dtype=int)
    coordinates = np.zeros((n, 3))

    for i, atom in enumerate(obmol.atoms):
        atomicnums[i] = atom.atomicnum
        coordinates[i] = atom.coords

    return Molecule(atomicnums, coordinates)


class Molecule:
    def __init__(self, atomicnums, coordinates):

        atomicnums = np.asarray(atomicnums, dtype=int)
        coordinates = np.asarray(coordinates, dtype=float)

        self.natoms = len(atomicnums)

        assert atomicnums.shape == (self.natoms,)
        assert coordinates.shape == (self.natoms, 3)

        self.atomicnums = atomicnums
        self.coordinates = coordinates

        self.stripped = np.all(atomicnums != 1)

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
