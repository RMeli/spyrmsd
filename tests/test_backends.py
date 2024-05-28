import numpy as np
import pytest

import spyrmsd
from spyrmsd import graph

# TODO: Run even with two backends installed


@pytest.mark.filterwarnings(
    "ignore::UserWarning"
)  # Silence "The backend is already" warning
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


def test_set_backend_unknown():
    with pytest.raises(ValueError, match="backend is not recognized or supported"):
        spyrmsd.set_backend("unknown")


def test_set_backend_same():
    current_backend = spyrmsd.get_backend()
    with pytest.warns(UserWarning, match=f"The backend is already {current_backend}."):
        spyrmsd.set_backend(current_backend)


@pytest.mark.filterwarnings(
    "ignore::UserWarning"
)  # Silence "The backend is already" warning
@pytest.mark.skipif(
    # Run test if all supported backends are installed
    not set(spyrmsd.graph._supported_backends) <= set(spyrmsd.available_backends),
    reason="Not all of the required backends are installed",
)
def test_molecule_graph_cache(mol) -> None:
    import graph_tool as gt
    import networkx as nx
    import rustworkx as rx

    m = mol.mol

    # Check molecules is in a clean state
    assert len(m.G.items()) == 0
    assert not m.stripped

    spyrmsd.set_backend("networkx")
    m.to_graph()

    assert "networkx" in m.G.keys()
    assert "graph_tool" not in m.G.keys()
    assert "rustworkx" not in m.G.keys()

    spyrmsd.set_backend("graph-tool")
    m.to_graph()

    assert "networkx" in m.G.keys()
    assert "graph_tool" in m.G.keys()
    assert "rustworkx" not in m.G.keys()

    spyrmsd.set_backend("rustworkx")
    m.to_graph()

    assert "networkx" in m.G.keys()
    assert "graph_tool" in m.G.keys()
    assert "rustworkx" in m.G.keys()

    # Make sure all backends (still) have a cache
    assert isinstance(m.G["networkx"], nx.Graph)
    assert isinstance(m.G["graph_tool"], gt.Graph)
    assert isinstance(m.G["rustworkx"], rx.PyGraph)

    # Strip molecule to ensure the cache is reset
    m.strip()

    assert len(m.G.items()) == 0
