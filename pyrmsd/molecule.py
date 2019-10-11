from pyrmsd import utils

from openbabel import pybel
import qcelemental as qcel

import numpy as np


def load(fname: str):

    fmt = utils.format_openbabel(fname)

    mol = next(pybel.readfile(fmt, fname))

    return mol


def openbabel_to_molecule(mol):

    n = len(mol.atoms)

    atomicnums = np.zeros(n, dtype=int)
    coordinates = np.zeros((n, 3))

    for i, atom in enumerate(mol.atoms):
        atomicnums[i] = atom.atomicnum
        coordinates[i] = atom.coords

    return Molecule(atomicnums, coordinates)


class Molecule:
    def __init__(self, atomicnums, coordinates):

        self.natoms = len(atomicnums)

        assert coordinates.shape == (self.natoms, 3)

        self.atomicnums = atomicnums
        self.coordinates = coordinates

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

    def to_graph(self):
        raise NotImplementedError

    def __len__(self):
        return self.natoms
