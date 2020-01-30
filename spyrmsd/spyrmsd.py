"""
Python RMSD tool
"""

import numpy as np

from spyrmsd import molecule, rmsd


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


def rmsdwrapper(
    mol1, mol2, symmetry=True, center=False, minimize=False, strip=False
) -> float:
    """
    Compute RMSD between two molecule.

    Parameters
    ----------
    mol1: molecule.Molecule
        Molecule 1
    mol2: molecule.Molecule
        Molecule 2
    symmetry: bool, optional
        Symmetry-corrected RMSD (using graph isomorphism)
    center: bool, optional
        Center molecules at origin
    minimize: bool, optional
        Minimised RMSD (using the quaternion polynomial method)
    strip: bool, optional
        Strip hydrogen atoms

    Returns
    -------
    float
        RMSD
    """

    if strip:
        mol1.strip()  # Does nothing if already stripped
        mol2.strip()

    if minimize:
        center = True

    c1 = coords_from_molecule(mol1, center)
    c2 = coords_from_molecule(mol2, center)

    if c1.shape != c2.shape:
        # TODO: Create specific exception
        raise ValueError("Molecules have different sizes.")

    RMSD = np.inf

    if minimize and symmetry:
        RMSD = rmsd.rmsd_qcp_isomorphic(
            c1, c2, mol1.adjacency_matrix, mol2.adjacency_matrix
        )
    elif minimize and not symmetry:
        RMSD = rmsd.rmsd_qcp(c1, c2, mol1.atomicnums, mol2.atomicnums)
    elif not minimize and symmetry:
        RMSD = rmsd.rmsd_isomorphic(
            c1, c2, mol1.adjacency_matrix, mol2.adjacency_matrix
        )
    elif not minimize and not symmetry:
        RMSD = rmsd.rmsd_standard(c1, c2, mol1.atomicnums, mol2.atomicnums)

    return RMSD


if __name__ == "__main__":

    from pyrmsd import io

    import argparse as ap
    import os

    parser = ap.ArgumentParser(description="Python RMSD tool.")

    parser.add_argument("reference", type=str, help="Reference file")
    parser.add_argument("molecules", type=str, nargs="+", help="Input file(s)")
    parser.add_argument("-m", "--minimize", action="store_true", help="Minimize (fit)")
    parser.add_argument(
        "-c", "--center", action="store_true", help="Center molecules at origin"
    )
    parser.add_argument("-s", "--strip", action="store_true", help="Strip H atoms")
    parser.add_argument(
        "-n", "--nosymm", action="store_false", help="No graph isomorphism"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    output: str = ""

    obref = io.load(args.reference)
    ref = io.to_molecule(obref, adjacency=True)

    if args.verbose:
        refname = os.path.basename(args.reference)

    # Loop over input files
    for molfile in args.molecules:

        # Load all molecule within file
        obmols = io.loadall(molfile)
        mols = [io.to_molecule(obmol, adjacency=True) for obmol in obmols]

        if args.verbose:
            molname = os.path.basename(molfile)

            output = f"{refname}:{molname} "

        # Loop over molecules within file
        for idx, mol in enumerate(mols):

            r = rmsdwrapper(
                ref, mol, args.nosymm, args.center, args.minimize, args.strip
            )

            print(f"{output}{r:.5f}")
