import numpy as np
import pytest

import spyrmsd
from spyrmsd import graph
from tests import molecules


@pytest.mark.skipif(
    # Run test if all supported backends are installed
    not set(spyrmsd.graph._supported_backends) <= set(spyrmsd.available_backends),
    reason="Not all of the required backends are installed",
)
def test_set_backend() -> None:
    import graph_tool as gt
    import networkx as nx
    import rustworkx as rx

    A = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 1]])

    spyrmsd.set_backend("networkx")
    assert spyrmsd.get_backend() == "networkx"

    Gnx = graph.graph_from_adjacency_matrix(A)
    assert isinstance(Gnx, nx.Graph)

    spyrmsd.set_backend("graph-tool")
    assert spyrmsd.get_backend() == "graph_tool"

    Ggt = graph.graph_from_adjacency_matrix(A)
    assert isinstance(Ggt, gt.Graph)

    spyrmsd.set_backend("rustworkx")
    assert spyrmsd.get_backend() == "rustworkx"

    Grx = graph.graph_from_adjacency_matrix(A)
    assert isinstance(Grx, rx.PyGraph)

    with pytest.raises(ValueError, match="backend is not recognized or supported"):
        spyrmsd.set_backend("unknown")


@pytest.mark.skipif(
    # Run test if all supported backends are installed
    not set(spyrmsd.graph._supported_backends) <= set(spyrmsd.available_backends),
    reason="Not all of the required backends are installed",
)
@pytest.mark.parametrize(
    "mol",
    [(molecules.benzene), (molecules.ethanol), (molecules.dialanine)],
    ids=["benzene", "ethanol", "dialanine"],
)
def test_molecule_graph_cache(mol) -> None:
    import graph_tool as gt
    import networkx as nx
    import rustworkx as rx

    ## Graph cache persists from previous tests, manually reset them
    mol.G = {}
    spyrmsd.set_backend("networkx")
    mol.to_graph()

    assert isinstance(mol.G["networkx"], nx.Graph)
    assert "graph_tool" not in mol.G.keys()

    spyrmsd.set_backend("graph-tool")
    mol.to_graph()

    spyrmsd.set_backend("rustworkx")
    mol.to_graph()

    ## Make sure all backends (still) have a cache
    assert isinstance(mol.G["networkx"], nx.Graph)
    assert isinstance(mol.G["graph_tool"], gt.Graph)
    assert isinstance(mol.G["rustworkx"], rx.PyGraph)

    ## Strip the molecule to ensure the cache is reset
    mol.strip()

    assert len(mol.G.items()) == 0
