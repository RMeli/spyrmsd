import networkx as nx
import qcelemental as qcel
import numpy as np

from openbabel import openbabel as ob

covalent_bond_multiplier: float = 1.2


def adjacency_matrix_from_obmol(obmol):

    n = len(obmol.atoms)

    A = np.zeros((n, n), dtype=int)

    for bond in ob.OBMolBondIter(obmol.OBMol):
        i = bond.GetBeginAtomIdx() - 1
        j = bond.GetEndAtomIdx() - 1

        A[i, j] = A[j, i] = 1

    return A


def graph_from_molecule(atomicnums, coordinates, named=False):

    n = len(atomicnums)

    assert coordinates.shape == (n, 3)

    G = nx.Graph()

    for i in range(n):
        r_i = qcel.covalentradii.get(atomicnums[i], units="angstrom")

        for j in range(i + 1, n):
            r_j = qcel.covalentradii.get(atomicnums[i], units="angstrom")

            distance = np.sqrt(np.sum((coordinates[i] - coordinates[j]) ** 2))

            if distance < (r_i + r_j) * covalent_bond_multiplier:
                G.add_edge(i, j)

    assert G.number_of_nodes() == n

    if named:
        mapping = {
            i: {"element": qcel.periodictable.to_symbol(anum)}
            for i, anum in enumerate(atomicnums)
        }
        nx.set_node_attributes(G, mapping)

    return G


if __name__ == "__main__":

    from pyrmsd import molecule

    import argparse as ap
    from matplotlib import pyplot as plt

    parser = ap.ArgumentParser(description="Draw molecular graph.")

    parser.add_argument("input", type=str, help="Input file")
    parser.add_argument("-o", "--output", type=str, default=None, help="Input file")

    args = parser.parse_args()

    mol = molecule.load(args.input)
    m = molecule.openbabel_to_molecule(mol)

    G = graph_from_molecule(m.atomicnums, m.coordinates, named=True)
    labels = nx.get_node_attributes(G, "element")

    nx.draw_kamada_kawai(G, labels=labels)
    if args.output is None:
        plt.plot()
    else:
        plt.savefig(args.output)
