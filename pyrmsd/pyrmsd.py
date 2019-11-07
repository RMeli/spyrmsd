"""
Python RMSD tool
"""

if __name__ == "__main__":

    from pyrmsd import molecule, rmsd

    import argparse as ap
    import os

    parser = ap.ArgumentParser(description="Python RMSD tool.")

    parser.add_argument("reference", type=str, help="Reference file")
    parser.add_argument("molecules", type=str, nargs="+", help="Input file")
    parser.add_argument("-s", "--strip", action="store_true", help="Strip H atoms")
    parser.add_argument("-m", "--minimize", action="store_true", help="Strip H atoms")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    output: str = ""

    obref = molecule.load(args.reference)
    ref = molecule.openbabel_to_molecule(obref, adjacency=True)

    if args.verbose:
        refname = os.path.basename(args.reference)

    # Loop over input files
    for molfile in args.molecules:

        # Load all molecule within file
        obmols = molecule.loadall(molfile)
        mols = [
            molecule.openbabel_to_molecule(obmol, adjacency=True) for obmol in obmols
        ]

        if args.verbose:
            molname = os.path.basename(molfile)
            output = f"{refname}:{molname} "

        # Loop over molecules within file
        for idx, mol in enumerate(mols):

            if args.strip:
                ref.strip()  # Does nothing if already stripped
                mol.strip()

            if args.minimize:  # QCP method
                r = rmsd.rmsd_qcp(ref, mol)
            else:  # Exact RMSD using graph isomorphism
                r = rmsd.rmsd_isomorphic(ref, mol, center=False)

            print(f"{output}{r:.5f}")
