from typing import Dict, List, Optional, Tuple

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


def _rmsd_isomorphic_core(
    coords1: np.ndarray,
    coords2: np.ndarray,
    am1: np.ndarray,
    am2: np.ndarray,
    atomicnums1: np.ndarray = None,
    atomicnums2: np.ndarray = None,
    center: bool = False,
    minimize: bool = False,
    isomorphisms: Optional[List[Dict[int, int]]] = None,
) -> Tuple[float, List[Dict[int, int]]]:
    """
    Compute RMSD using graph isomorphism.

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
    atomicnums1: npndarray, optional
        Atomic numbers for molecule 1
    atomicnums2: npndarray, optional
        Atomic numbers for molecule 2
    center: bool
        Centering flag
    minimize: bool
        Compute minized RMSD
    isomorphisms: Optional[List[Dict[int,int]]]
        Previously computed graph isomorphism

    Returns
    -------
    Tuple[float, List[Dict[int, int]]]
        RMSD (after graph matching) and graph isomorphisms
    """

    assert coords1.shape == coords2.shape

    n = coords1.shape[0]

    # Center coordinates if required
    c1 = utils.center(coords1) if center or minimize else coords1
    c2 = utils.center(coords2) if center or minimize else coords2

    # No cached isomorphisms
    if isomorphisms is None:
        # Convert molecules to graphs
        G1 = graph.graph_from_adjacency_matrix(am1, atomicnums1)
        G2 = graph.graph_from_adjacency_matrix(am2, atomicnums2)

        # Get all the possible graph isomorphisms
        isomorphisms = graph.match_graphs(G1, G2)

    # Minimum result
    # Squared displacement (not minimize) or RMSD (minimize)
    min_result = np.inf

    # Loop over all graph isomorphisms to find the lowest RMSD
    for isomorphism in isomorphisms:

        # Use the isomorphism to shuffle coordinates around (from original order)
        c1i = c1[list(isomorphism.keys()), :]
        c2i = c2[list(isomorphism.values()), :]

        if not minimize:
            # Compute square displacement
            # Avoid dividing by n and an expensive sqrt() operation
            result = np.sum((c1i - c2i) ** 2)
        else:
            # Compute minimized RMSD using QCP
            result = qcp.qcp_rmsd(c1i, c2i)

        min_result = result if result < min_result else min_result

    if not minimize:
        # Compute actual RMSD from square displacement
        min_result = np.sqrt(min_result / n)

    # Return the actual RMSD
    return min_result, isomorphisms


def rmsd_isomorphic(
    coords1: np.ndarray,
    coords2: np.ndarray,
    am1: np.ndarray,
    am2: np.ndarray,
    atomicnums1: np.ndarray = None,
    atomicnums2: np.ndarray = None,
    center: bool = False,
    minimize: bool = False,
) -> float:
    """
    Compute RMSD using graph isomorphism.

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
    atomicnums1: npndarray, optional
        Atomic numbers for molecule 1
    atomicnums2: npndarray, optional
        Atomic numbers for molecule 2
    center: bool
        Centering flag
    minimize: bool
        Minimum RMSD

    Returns
    -------
    float
        RMSD (after graph matching) and graph isomorphisms

    Notes
    -----

    This QCP method, activated with the keyword `minimize=True` works in cases where
    the atoms in `mol1` and `mol2` are not in the exact same order. If the atoms in
    `mol1` and `mol2` are in the same order use `rmsd_qcp` (faster).
    """

    RMSD, _ = _rmsd_isomorphic_core(
        coords1,
        coords2,
        am1,
        am2,
        atomicnums1,
        atomicnums2,
        center=center,
        minimize=minimize,
        isomorphisms=None,
    )

    return RMSD


def multirmsd_isomorphic(
    coordsref: np.ndarray,
    coords: List[np.ndarray],
    amref: np.ndarray,
    am: np.ndarray,
    atomicnumsref: np.ndarray = None,
    atomicnums: np.ndarray = None,
    center: bool = False,
    minimize: bool = False,
    cache: bool = True,
) -> List[float]:
    """
    Compute RMSD using graph isomorphism for multiple coordinates.

    Parameters
    ----------
    coordsref: np.ndarray
        Coordinate of reference molecule
    coords: List[np.ndarray]
        Coordinates of other molecule
    amref: np.ndarray
        Adjacency matrix for reference molecule
    am: np.ndarray
        Adjacency matrix for other molecule
    atomicnumsref: npndarray, optional
        Atomic numbers for reference
    atomicnums: npndarray, optional
        Atomic numbers for other molecule
    center: bool
        Centering flag
    minimize: bool
        Minimum RMSD

    Returns
    -------
    float
        RMSD (after graph matching) and graph isomorphisms

    Notes
    -----

    This QCP method, activated with the keyword `minimize=True` works in cases where
    the atoms in `mol1` and `mol2` are not in the exact same order. If the atoms in
    `mol1` and `mol2` are in the same order use `rmsd_qcp` (faster).
    """

    RMSDlist, isomorphism = [], None

    for c in coords:

        if not cache:
            # Reset isomorphism
            isomorphism = None

        RMSD, isomorphism = _rmsd_isomorphic_core(
            coordsref,
            c,
            amref,
            am,
            atomicnumsref,
            atomicnums,
            center=center,
            minimize=minimize,
            isomorphisms=isomorphism,
        )

        RMSDlist.append(RMSD)

    return RMSDlist
