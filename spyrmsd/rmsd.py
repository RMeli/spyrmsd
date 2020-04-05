from typing import Any, List, Optional, Tuple, Union

import numpy as np

from spyrmsd import graph, hungarian, qcp, utils


def rmsd(
    coords1: np.ndarray,
    coords2: np.ndarray,
    atomicn1: np.ndarray,
    atomicn2: np.ndarray,
    center: bool = False,
    minimize: bool = False,
) -> float:
    """
    Compute RMSD

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
    center: bool
        Center molecules at origin
    minimize: bool
        Compute minimum RMSD (with QCP method)

    Returns
    -------
    float
        RMSD

    Notes
    -----
    When `minimize=True`, the QCP method is used. [1]_ The molecules are
    centred at the origin according to the center of geometry and superimposed
    in order to minimize the RMSD.

    .. [1] D. L. Theobald, *Rapid calculation of RMSDs using a quaternion-based
       characteristic polynomial*, Acta Crys. A **61**, 478-480 (2005).
    """

    assert np.all(atomicn1 == atomicn2)
    assert coords1.shape == coords2.shape

    # Center coordinates if required
    c1 = utils.center(coords1) if center or minimize else coords1
    c2 = utils.center(coords2) if center or minimize else coords2

    if minimize:
        rmsd = qcp.qcp_rmsd(c1, c2)
    else:
        n = coords1.shape[0]

        rmsd = np.sqrt(np.sum((c1 - c2) ** 2) / n)

    return rmsd


def hrmsd(
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
    atomicnums1: np.ndarray,
    atomicnums2: np.ndarray,
    am1: np.ndarray,
    am2: np.ndarray,
    center: bool = False,
    minimize: bool = False,
    isomorphisms: Optional[List[Tuple[List[int], List[int]]]] = None,
) -> Tuple[float, List[Tuple[List[int], List[int]]]]:
    """
    Compute RMSD using graph isomorphism.

    Parameters
    ----------
    coords1: np.ndarray
        Coordinate of molecule 1
    coords2: np.ndarray
        Coordinates of molecule 2
    atomicnums1: npndarray
        Atomic numbers for molecule 1
    atomicnums2: npndarray
        Atomic numbers for molecule 2
    am1: np.ndarray
        Adjacency matrix for molecule 1
    am2: np.ndarray
        Adjacency matrix for molecule 2
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
    for idx1, idx2 in isomorphisms:

        # Use the isomorphism to shuffle coordinates around (from original order)
        c1i = c1[idx1, :]
        c2i = c2[idx2, :]

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


def symmrmsd(
    coordsref: np.ndarray,
    coords: Union[np.ndarray, List[np.ndarray]],
    atomicnumsref: np.ndarray,
    atomicnums: np.ndarray,
    amref: np.ndarray,
    am: np.ndarray,
    center: bool = False,
    minimize: bool = False,
    cache: bool = True,
) -> Any:
    """
    Compute RMSD using graph isomorphism for multiple coordinates.

    Parameters
    ----------
    coordsref: np.ndarray
        Coordinate of reference molecule
    coords: List[np.ndarray]
        Coordinates of other molecule
    atomicnumsref: npndarray
        Atomic numbers for reference
    atomicnums: npndarray
        Atomic numbers for other molecule
    amref: np.ndarray
        Adjacency matrix for reference molecule
    am: np.ndarray
        Adjacency matrix for other molecule
    center: bool
        Centering flag
    minimize: bool
        Minimum RMSD

    Returns
    -------
    float: Union[float, List[float]]
        Symmetry-corrected RMSD(s) and graph isomorphisms

    Notes
    -----

    Graph isomorphism is introduced for symmetry corrections. However, it is also
    useful when two molecules do not have the atoms in the same order since atom
    matching according to atomic numbers and the molecular connectivity is
    performed. If atoms are in the same order and there is no symmetry, use the
    `rmsd` function.
    """

    if isinstance(coords, list):  # Multiple RMSD calculations

        RMSD: Any = []
        isomorphism = None

        for c in coords:

            if not cache:
                # Reset isomorphism
                isomorphism = None

            srmsd, isomorphism = _rmsd_isomorphic_core(
                coordsref,
                c,
                atomicnumsref,
                atomicnums,
                amref,
                am,
                center=center,
                minimize=minimize,
                isomorphisms=isomorphism,
            )

            RMSD.append(srmsd)

    else:  # Single RMSD calculation

        RMSD, isomorphism = _rmsd_isomorphic_core(
            coordsref,
            coords,
            atomicnumsref,
            atomicnums,
            amref,
            am,
            center=center,
            minimize=minimize,
            isomorphisms=None,
        )

    return RMSD
