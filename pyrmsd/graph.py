import networkx as nx
import qcelemental as qcel
import numpy as np

covalent_bond_multiplier: float = 1.3

def graph_from_molecule(atomicnums, coordinates):

    n = len(atomicnums)

    assert coordinates.shape == (n, 3)

    G = nx.Graph()

    for i in range(n):
        r_i = qcel.covalentradii.get(atomicnums[i], units='angstrom')

        for j in range(i + 1, n):
            r_j = qcel.covalentradii.get(atomicnums[i], units='angstrom')

            distance = np.sqrt(np.sum((coordinates[i] - coordinates[j])**2))

            if distance < (r_i + r_j) * covalent_bond_multiplier:
                G.add_edge(i, j)

    assert G.number_of_nodes() == n

    return G

if __name__ == "__main__":

    from pyrmsd import molecule, graph

    import argparse as ap
    from matplotlib import pyplot as plt

    parser = ap.ArgumentParser(description="Draw molecular graph.")

    parser.add_argument("input", type=str, help="Input file")
    parser.add_argument("-o", "--output", type=str, default=None, help="Input file")

    args = parser.parse_args()

    mol = molecule.load(args.input)
    m = molecule.openbabel_to_molecule(mol)

    G = graph_from_molecule(m.atomicnums, m.coordinates)

    nx.draw_kamada_kawai(G)
    if args.output is None:
        plt.plot()
    else:
        plt.savefig(args.output)