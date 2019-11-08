from pyrmsd import qcp, hungarian, graph, molecule

import numpy as np


def coords_from_molecule(mol: molecule.Molecule, center: bool = False) -> np.ndarray:
    """
    Atomic coordinates from molecule.

    Parameters
    ----------
    mol: molecule.Molecule
        Molecule
    center: bool
        Center flag

    Returns
    -------
    np.ndarray
        Atomic coordinates (possibly centred)

    Notes
    -----
    Atomic coordinates are centred according to the center of geometry, not the center
    of mass.
    """

    if center:
        coords = mol.coordinates - mol.center_of_geometry()
    else:
        coords = mol.coordinates

    return coords


def rmsd_dummy(
    mol1: molecule.Molecule, mol2: molecule.Molecule, center: bool = False
) -> float:
    """
    Compute dummy (naïve) RMSD.

    Parameters
    ----------
    mol1: molecule.Molecule
        Molecule 1
    mol2: molecule.Molecule
        Molecule 2
    center: bool
        Flag for centering the molecules at the origin

    Returns
    -------
    float
        Dummy (naïve) RMSD

    Notes
    -----
    The dummy (naïve) RMSD distribution assumes a one-to-one mapping of the atoms in
    the order they are stored.

    The equality between all the atomic numbers is checked for sanity.
    """

    assert np.all(mol1.atomicnums == mol2.atomicnums)
    assert mol1.coordinates.shape == mol2.coordinates.shape

    n = len(mol1)

    c1 = coords_from_molecule(mol1, center)
    c2 = coords_from_molecule(mol2, center)

    return np.sqrt(np.sum((c1 - c2) ** 2) / n)


def rmsd_qcp(mol1, mol2) -> float:
    """
    Compute minimum RMSD using the Quaternion Characteristic Polynomial method.

    Parameters
    ----------
    mol1: molecule.Molecule
        Molecule 1
    mol2: molecule.Molecule
        Molecule 2

    Returns
    -------
    float
        Minimum RMSD (after superimposition)

    Notes
    -----
    The molecules are always centred at the origin according to the center of geometry
    and superimposed in order to minimize the RMSD. See [1]_ for details.

    .. [1] D. L. Theobald, *Rapid calculation of RMSDs using a quaternion-based
       characteristic polynomial*, Acta Crys. A**61**, 478-480 (2005).
    """

    assert np.all(mol1.atomicnums == mol2.atomicnums)

    c1 = coords_from_molecule(mol1, center=True)
    c2 = coords_from_molecule(mol2, center=True)

    return qcp.qcp_rmsd(c1, c2)


def rmsd_hungarian(mol1, mol2, center=False):
    """
    Compute minimum RMSD using the Hungarian method.

    Parameters
    ----------
    mol1: molecule.Molecule
        Molecule 1
    mol2: molecule.Molecule
        Molecule 2

    Returns
    -------
    float
        Minimum RMSD (after assignment)

    Notes
    -----
    The Hungarian algorithm is used to solve the linear assignment problem, which is
    a minimum weight matching of the molecular graphs (bipartite). See [2]_ for details.

    The linear assignment problem is solved for every element separately.

    .. [2] W. J. Allen and R. C. Rizzo, *Implementation of the Hungarian Algorithm to
        Account for Ligand Symmetry and Similarity in Structure-Based Design*,
        J. Chem. Inf. Model. **54**, 518-529 (2014)
    """

    assert mol1.atomicnums.shape == mol2.atomicnums.shape
    assert mol1.coordinates.shape == mol2.coordinates.shape

    c1 = coords_from_molecule(mol1, center)
    c2 = coords_from_molecule(mol2, center)

    return hungarian.hungarian_rmsd(c1, c2, mol1.atomicnums, mol2.atomicnums)


def rmsd_isomorphic(mol1, mol2, center=False):
    """
    Compute minimum RMSD using graph isomorphism.

    Parameters
    ----------
    mol1: molecule.Molecule
        Molecule 1
    mol2: molecule.Molecule
        Molecule 2

    Returns
    -------
    float
        Minimum RMSD (after graph matching)
    """

    assert mol1.atomicnums.shape == mol2.atomicnums.shape
    assert mol1.coordinates.shape == mol2.coordinates.shape

    n = len(mol1)

    c1 = coords_from_molecule(mol1, center)
    c2 = coords_from_molecule(mol2, center)

    # Convert molecules to graphs
    G1 = mol1.to_graph()
    G2 = mol2.to_graph()

    # Get all the possible graph isomorphisms
    isomorphisms = graph.match_graphs(G1, G2)

    # Minimum squared displacement
    min_sd = np.inf

    # Loop over all graph isomorphisms to find the lowest RMSD
    for isomorphism in isomorphisms:

        # Use the isomorphism to shuffle coordinates around (from original order)
        c1i = c1[list(isomorphism.keys()), :]
        c2i = c2[list(isomorphism.values()), :]

        # Compute square displacement
        # Avoid dividing by n and an expensive sqrt() operation
        sd = np.sum((c1i - c2i) ** 2)

        if sd < min_sd:
            min_sd = sd

    # Return the actual RMSD
    return np.sqrt(min_sd / n)


def rmsd_qcp_isomorphic(mol1: molecule.Molecule, mol2: molecule.Molecule) -> float:
    """
    Compute minimum RMSD using the Quaternion Characteristic Polynomial method on
    isomorphic graphs.

    Parameters
    ----------
    mol1: molecule.Molecule
        Molecule 1
    mol2: molecule.Molecule
        Molecule 2

    Returns
    -------
    float
        Minimum RMSD (after graph matching and superimposition)

    Notes
    -----
    This QCP method works in cases where the atoms in `mol1` and `mol2` are not in the
    exact same order. If the atoms in `mol1` and `mol2` are in the same order use
    `rmsd_qcp` (faster).
    """

    assert mol1.atomicnums.shape == mol2.atomicnums.shape
    assert mol1.coordinates.shape == mol2.coordinates.shape

    c1 = coords_from_molecule(mol1, center=True)
    c2 = coords_from_molecule(mol2, center=True)

    # Convert molecules to graphs
    G1 = mol1.to_graph()
    G2 = mol2.to_graph()

    # Get all the possible graph isomorphisms
    isomorphisms = graph.match_graphs(G1, G2)

    # Minimum squared displacement
    min_rmsd = np.inf

    # Loop over all graph isomorphisms to find the lowest RMSD
    for isomorphism in isomorphisms:

        # Use the isomorphism to shuffle coordinates around (from original order)
        c1i = c1[list(isomorphism.keys()), :]
        c2i = c2[list(isomorphism.values()), :]

        rmsd = qcp.qcp_rmsd(c1i, c2i)

        if rmsd < min_rmsd:
            min_rmsd = rmsd

    return min_rmsd
