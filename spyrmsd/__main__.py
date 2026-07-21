"""
Symmetry-corrected RMSD calculations in Python
"""

import argparse as ap
import importlib.util
import sys
import warnings

import spyrmsd
from spyrmsd import graph, io, molecule
from spyrmsd.rmsd import rmsdwrapper


def build_parser() -> ap.ArgumentParser:
    parser = ap.ArgumentParser(
        prog="spyrmsd",
        description="Symmetry-corrected RMSD calculations in Python.",
    )

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
    parser.add_argument(
        "--graph-builder",
        dest="graph_builder",
        choices=("simple", "xyzgraph"),
        default="simple",
        help=(
            "Adjacency matrix builder "
            "(xyzgraph: Molecular Graph Construction from Cartesian Coordinates)"
        ),
    )
    parser.add_argument(
        "-g",
        "--graph-backend",
        type=str,
        default=None,
        help="Graph library (backend)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose mode"
    )
    parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {spyrmsd.__version__}"
    )

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if importlib.util.find_spec("rdkit") is None:
        raise ImportError(
            "RDKit not found. Please install RDKit to use sPyRMSD as a standalone tool."
        )

    try:
        ref = io.loadmol(args.reference)
    except OSError:
        print("ERROR: Reference file not found.", file=sys.stderr)
        return -1

    try:
        mols = [mol for molfile in args.molecules for mol in io.loadallmols(molfile)]
    except OSError:
        print("ERROR: Molecule file(s) not found.", file=sys.stderr)
        return -1

    if args.graph_backend is not None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            spyrmsd.set_backend(args.graph_backend)

    try:
        graph.set_graph_builder(args.graph_builder)
        ref = molecule.Molecule(
            ref.atomicnums,
            ref.coordinates,
            adjacency_matrix=graph.adjacency_matrix_from_atomic_coordinates(
                ref.atomicnums, ref.coordinates
            ),
        )
        mols = [
            molecule.Molecule(
                mol.atomicnums,
                mol.coordinates,
                adjacency_matrix=graph.adjacency_matrix_from_atomic_coordinates(
                    mol.atomicnums, mol.coordinates
                ),
            )
            for mol in mols
        ]
    except (ImportError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return -1

    if args.verbose:
        print(f"Graph library: {spyrmsd.get_backend()}")
        print(f"Graph builder: {graph.get_graph_builder()}")

    rmsd_list = rmsdwrapper(
        ref,
        mols,
        symmetry=args.nosymm,
        center=args.center,
        minimize=args.minimize,
        strip=not args.hydrogens,
    )

    for rmsd in rmsd_list:
        print(f"{rmsd:.5f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
