from pyrmsd import utils, graph

import qcelemental as qcel
import numpy as np

import networkx as nx

from typing import List


class Molecule:
    def __init__(self, atomicnums, coordinates, adjacency_matrix=None):

        atomicnums = np.asarray(atomicnums, dtype=int)
        coordinates = np.asarray(coordinates, dtype=float)

        self.natoms: int = len(atomicnums)

        assert atomicnums.shape == (self.natoms,)
        assert coordinates.shape == (self.natoms, 3)

        self.atomicnums = atomicnums
        self.coordinates = coordinates

        self.stripped: bool = np.all(atomicnums != 1)

        if adjacency_matrix is not None:
            self.adjacency_matrix = np.asarray(adjacency_matrix, dtype=int)

        self.G: nx.Graph = None

        self.masses: List[float] = None

    def translate(self, vector: np.ndarray):
        """
        Translate molecule.

        Parameters
        ----------
        vector: np.ndarray
            Translation vector (in 3D)
        """
        assert len(vector) == 3
        self.coordinates += vector

    def rotate(self, angle: float, axis: np.ndarray, units: str = "rad"):
        """
        Rotate molecule.

        Parameters
        ----------
        angle: float
            Rotation angle
        axis: np.ndarray
            Axis of rotation (in 3D)
        units: {"rad", "deg"}
            Units of the angle (radians `rad` or degrees `deg`)
        """
        assert len(axis) == 3
        self.coordinates = utils.rotate(self.coordinates, angle, axis, units)

    def center_of_mass(self) -> np.ndarray:
        """
        Center of mass

        Returns
        -------
        np.ndarray
            Center of mass

        Notes
        -----
        Atomic masses are cached.
        """
        if self.masses is None:
            self.masses = [qcel.periodictable.to_mass(anum) for anum in self.atomicnums]

        return np.average(self.coordinates, axis=0, weights=self.masses)

    def center_of_geometry(self) -> np.ndarray:
        """
        Center of geometry

        Returns
        -------
        np.ndarray
            Center of geometry
        """
        return utils.center_of_geometry(self.coordinates)

    def strip(self) -> None:
        """
        Strip hydrogen atoms.
        """

        if not self.stripped:
            idx = self.atomicnums != 1  # Hydrogen atoms

            # Strip
            self.atomicnums = self.atomicnums[idx]
            self.coordinates = self.coordinates[idx, :]

            # Update number of atoms
            self.natoms = len(self.atomicnums)

            # Update adjacency matrix
            self.adjacency_matrix = self.adjacency_matrix[np.ix_(idx, idx)]

            self.stripped = True

    def to_graph(self) -> nx.Graph:
        """
        Convert molecule to graph.

        Returns
        -------
        networkx.Graph
            Molecular graph.

        Raises
        ------
        NotImplementedError
            If there is no associated adjacency matrix.

        Notes
        -----
        The molecule needs to have an associated adjacency matrix. Bond perception
        will be implemented in later versions.

        The molecular graph is cached.
        """
        if self.G is None:
            if self.adjacency_matrix is not None:
                self.G = graph.graph_from_adjacency_matrix(self.adjacency_matrix)
            else:
                raise NotImplementedError

        return self.G

    def __len__(self) -> int:
        """
        Molecule size.

        Returns
        -------
        int
            Number of atoms within the molecule
        """
        return self.natoms
