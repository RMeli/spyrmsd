from pyrmsd import molecule, graph
from pyrmsd.tests import molecules


def test_graph_from_molecule_benzene_xyz():

    mol = molecules.benzene

    m = molecule.openbabel_to_molecule(mol)

    G = graph.graph_from_molecule(m.atomicnums, m.coordinates)

    assert G.number_of_nodes() == len(m)
    assert G.number_of_edges() == 12


def test_graph_from_molecule_ethanol_xyz():

    mol = molecules.ethanol

    m = molecule.openbabel_to_molecule(mol)

    G = graph.graph_from_molecule(m.atomicnums, m.coordinates)

    assert G.number_of_nodes() == len(m)
    assert G.number_of_edges() == 8


def test_graph_from_molecule_dialanine_sdf():

    mol = molecules.dialanine

    m = molecule.openbabel_to_molecule(mol)

    G = graph.graph_from_molecule(m.atomicnums, m.coordinates)

    assert G.number_of_nodes() == len(m)
    assert G.number_of_edges() == 22
