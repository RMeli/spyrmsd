from pyrmsd import graph
from pyrmsd.tests import molecules


def test_graph_from_molecule_benzene() -> None:

    mol = molecules.benzene

    G = graph.graph_from_molecule(mol.atomicnums, mol.coordinates)

    assert G.number_of_nodes() == len(mol)
    assert G.number_of_edges() == 12


def test_graph_from_molecule_named_benzene() -> None:

    mol = molecules.benzene

    G = graph.graph_from_molecule(mol.atomicnums, mol.coordinates, named=True)

    assert G.number_of_nodes() == len(mol)
    assert G.number_of_edges() == 12

    nodes = G.nodes()
    for node in nodes:
        element = nodes[node]["element"]

        assert element == "H" or element == "C"


def test_graph_from_molecule_ethanol() -> None:

    mol = molecules.ethanol

    G = graph.graph_from_molecule(mol.atomicnums, mol.coordinates)

    assert G.number_of_nodes() == len(mol)
    assert G.number_of_edges() == 8


def test_graph_from_molecule_dialanine() -> None:

    mol = molecules.dialanine

    G = graph.graph_from_molecule(mol.atomicnums, mol.coordinates)

    assert G.number_of_nodes() == len(mol)
    assert G.number_of_edges() == 22
