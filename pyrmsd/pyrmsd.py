"""
Python RMSD tool
"""

from pyrmsd import molecule, rmsd

if __name__ == "__main__":

    import argparse as ap

    parser = ap.ArgumentParser(description="Python RMSD tool.")

    parser.add_argument("reference", type=str, help="Reference file")
    parser.add_argument("molecule", type=str, help="Input file")
    parser.add_argument("-s", "--strip", action="store_true", help="Strip H atoms")
    parser.add_argument("-m", "--minimize", action="store_true", help="Strip H atoms")

    args = parser.parse_args()

    obref = molecule.load(args.reference)
    obmols = molecule.loadall(args.molecule)

    ref = molecule.openbabel_to_molecule(obref, adjacency=True)
    mols = [molecule.openbabel_to_molecule(obmol, adjacency=True) for obmol in obmols]

    for mol in mols:

        if args.strip:
            ref.strip()
            mol.strip()

        if args.minimize:
            r = rmsd.rmsd_qcp(ref, mol)
        else:
            r = rmsd.rmsd_isomorphic(ref, mol, center=False)

        print(f"{r:.5f}")
