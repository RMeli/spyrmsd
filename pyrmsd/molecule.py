from pyrmsd import utils

from openbabel import pybel
import qcelemental as qcel

import numpy as np

def load(fname: str):

    fmt = utils.format_openbabel(fname)

    mol = next(pybel.readfile(fmt, fname))

    return  mol

def openbabel_to_molecule(mol):

    n = len(mol.atoms)

    atomicnums = np.zeros(n, dtype=int)
    coordinates = np.zeros((n, 3))

    for i, atom in enumerate(mol.atoms):
        atomicnums[i] = atom.atomicnum
        coordinates[i] = atom.coords

    return Molecule(atomicnums, coordinates)
    
class Molecule():

    def __init__(self, atomicnums, coordinates):

        self.natoms = len(atomicnums)

        assert coordinates.shape == (self.natoms, 3)

        self.atomicnums = atomicnums
        self.coordinates = coordinates

    def to_graph():
        pass

    def __len__(self):
        return self.natoms


