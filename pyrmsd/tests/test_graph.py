from pyrmsd import molecule, graph
from pyrmsd.tests import molecules


def test_graph_from_molecule_benzene():

    mol = molecules.benzene

    m = molecule.openbabel_to_molecule(mol)

    G = graph.graph_from_molecule(m.atomicnums, m.coordinates)

    assert G.number_of_nodes() == len(m)
    assert G.number_of_edges() == 12


def test_graph_from_molecule_named_benzene():

    mol = molecules.benzene

    m = molecule.openbabel_to_molecule(mol)

    G = graph.graph_from_molecule(m.atomicnums, m.coordinates, named=True)

    assert G.number_of_nodes() == len(m)
    assert G.number_of_edges() == 12

    nodes = G.nodes()
    for node in nodes:
        element = nodes[node]["element"]

        assert element == "H" or element == "C"


def test_graph_from_molecule_ethanol():

    mol = molecules.ethanol

    m = molecule.openbabel_to_molecule(mol)

    G = graph.graph_from_molecule(m.atomicnums, m.coordinates)

    assert G.number_of_nodes() == len(m)
    assert G.number_of_edges() == 8


def test_graph_from_molecule_dialanine():

    mol = molecules.dialanine

    m = molecule.openbabel_to_molecule(mol)

    G = graph.graph_from_molecule(m.atomicnums, m.coordinates)

    assert G.number_of_nodes() == len(m)
    assert G.number_of_edges() == 22
