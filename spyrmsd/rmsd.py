import numpy as np

from spyrmsd import graph, hungarian, qcp, utils


def rmsd_standard(
    coords1: np.ndarray,
    coords2: np.ndarray,
    atomicn1: np.ndarray,
    atomicn2: np.ndarray,
    center: bool = False,
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

    assert np.all(atomicn1 == atomicn2)
    assert coords1.shape == coords2.shape

    n = coords1.shape[0]

    # Center coordinates if required
    c1 = utils.center(coords1) if center else coords1
    c2 = utils.center(coords2) if center else coords2

    return np.sqrt(np.sum((c1 - c2) ** 2) / n)


def rmsd_qcp(
    coords1: np.ndarray, coords2: np.ndarray, atomicn1: np.ndarray, atomicn2: np.ndarray
) -> float:
    """
    Compute minimum RMSD using the Quaternion Characteristic Polynomial method.

    Parameters
    ----------
    coords1: np.ndarray
        Coordinate of molecule 1
    coords2: np.ndarray
        Coordinates of molecule 2
    atomicn1: np.ndarray
        Atomic numbers for molecule 1
    atomicn2: np.ndarray
        Atomic numbers for molecule 2

    Returns
    -------
    float
        Minimum RMSD (after superimposition)

    Notes
    -----
    The molecules are always centred at the origin according to the center of geometry
    and superimposed in order to minimize the RMSD. [1]_

    .. [1] D. L. Theobald, *Rapid calculation of RMSDs using a quaternion-based
       characteristic polynomial*, Acta Crys. A **61**, 478-480 (2005).
    """

    assert np.all(atomicn1 == atomicn2)

    # Center coordinates if required
    c1 = utils.center(coords1)
    c2 = utils.center(coords2)

    return qcp.qcp_rmsd(c1, c2)


def rmsd_hungarian(
    coords1: np.ndarray,
    coords2: np.ndarray,
    atomicn1: np.ndarray,
    atomicn2: np.ndarray,
    center=False,
):
    """
    Compute minimum RMSD using the Hungarian method.

    Parameters
    ----------
    coords1: np.ndarray
        Coordinate of molecule 1
    coords2: np.ndarray
        Coordinates of molecule 2
    atomicn1: np.ndarray
        Atomic numbers for molecule 1
    atomicn2: np.ndarray
        Atomic numbers for molecule 2

    Returns
    -------
    float
        Minimum RMSD (after assignment)

    Notes
    -----
    The Hungarian algorithm is used to solve the linear assignment problem, which is
    a minimum weight matching of the molecular graphs (bipartite). [2]_

    The linear assignment problem is solved for every element separately.

    .. [2] W. J. Allen and R. C. Rizzo, *Implementation of the Hungarian Algorithm to
        Account for Ligand Symmetry and Similarity in Structure-Based Design*,
        J. Chem. Inf. Model. **54**, 518-529 (2014)
    """

    assert atomicn1.shape == atomicn2.shape
    assert coords1.shape == coords2.shape

    # Center coordinates if required
    c1 = utils.center(coords1) if center else coords1
    c2 = utils.center(coords2) if center else coords2

    return hungarian.hungarian_rmsd(c1, c2, atomicn1, atomicn2)


def rmsd_isomorphic(
    coords1: np.ndarray,
    coords2: np.ndarray,
    am1: np.ndarray,
    am2: np.ndarray,
    center=False,
) -> float:
    """
    Compute minimum RMSD using graph isomorphism.

    Parameters
    ----------
    coords1: np.ndarray
        Coordinate of molecule 1
    coords2: np.ndarray
        Coordinates of molecule 2
    am1: np.ndarray
        Adjacency matrix for molecule 1
    am2: np.ndarray
        Adjacency matrix for molecule 2
    center: boolean
        Centering flag

    Returns
    -------
    float
        Minimum RMSD (after graph matching)
    """

    assert coords1.shape == coords2.shape

    n = coords1.shape[0]

    # Center coordinates if required
    c1 = utils.center(coords1) if center else coords1
    c2 = utils.center(coords2) if center else coords2

    # Convert molecules to graphs
    G1 = graph.graph_from_adjacency_matrix(am1)
    G2 = graph.graph_from_adjacency_matrix(am2)

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


def rmsd_qcp_isomorphic(
    coords1: np.ndarray, coords2: np.ndarray, am1: np.ndarray, am2: np.ndarray
) -> float:
    """
    Compute minimum RMSD using the Quaternion Characteristic Polynomial method on
    isomorphic graphs.

    Parameters
    ----------
    coords1: np.ndarray
        Coordinate of molecule 1
    coords2: np.ndarray
        Coordinates of molecule 2
    am1: np.ndarray
        Adjacency matrix for molecule 1
    am2: np.ndarray
        Adjacency matrix for molecule 2

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

    assert coords1.shape == coords2.shape

    # Center coordinates
    c1 = utils.center(coords1)
    c2 = utils.center(coords2)

    # Build graph from adjacency matrix
    G1 = graph.graph_from_adjacency_matrix(am1)
    G2 = graph.graph_from_adjacency_matrix(am2)

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
