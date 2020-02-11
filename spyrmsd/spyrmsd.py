"""
Python RMSD tool
"""

from typing import List

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
    molref,
    mols,
    symmetry: bool = True,
    center: bool = False,
    minimize: bool = False,
    strip: bool = False,
    cache: bool = True,
) -> List[float]:
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
    List[float]
        RMSDs
    """

    if strip:
        molref.strip()

        for mol in mols:
            mol.strip()

    if minimize:
        center = True

    cref = coords_from_molecule(molref, center)
    cmols = [coords_from_molecule(mol, center) for mol in mols]

    RMSDlist = []

    if symmetry:
        RMSDlist = rmsd.multirmsd_isomorphic(
            cref,
            cmols,
            molref.adjacency_matrix,
            mols[0].adjacency_matrix,
            molref.atomicnums,
            mols[0].atomicnums,
            center=center,
            minimize=minimize,
            cache=cache,
        )
    elif minimize and not symmetry:
        for c in cmols:
            RMSDlist.append(
                rmsd.rmsd_qcp(cref, c, molref.atomicnums, mols[0].atomicnums)
            )
    elif not minimize and not symmetry:
        for c in cmols:
            RMSDlist.append(
                rmsd.rmsd_standard(cref, c, molref.atomicnums, mols[0].atomicnums)
            )

    return RMSDlist


if __name__ == "__main__":

    from spyrmsd import io

    import argparse as ap

    parser = ap.ArgumentParser(description="Python RMSD tool.")

    parser.add_argument("reference", type=str, help="Reference file")
    parser.add_argument("molecules", type=str, nargs="+", help="Input file(s)")
    parser.add_argument("-m", "--minimize", action="store_true", help="Minimize (fit)")
    parser.add_argument(
        "-c", "--center", action="store_true", help="Center molecules at origin"
    )
    parser.add_argument("--hydrogens", action="store_true", help="Keep hydrogen atoms")
    parser.add_argument(
        "-n", "--nosymm", action="store_false", help="No graph isomorphism"
    )

    args = parser.parse_args()

    inref = io.load(args.reference)
    ref = io.to_molecule(inref, adjacency=True)

    # Load all molecules
    inmols = [inmol for molfile in args.molecules for inmol in io.loadall(molfile)]
    mols = [io.to_molecule(mol, adjacency=True) for mol in inmols]

    # Loop over molecules within fil
    RMSDlist = rmsdwrapper(
        ref,
        mols,
        symmetry=args.nosymm,  # args.nosymm store False
        center=args.center,
        minimize=args.minimize,
        strip=not args.hydrogens,
    )

    for RMSD in RMSDlist:
        print(f"{RMSD:.5f}")
