import builtins
import sys
import types

import numpy as np
import pytest

from spyrmsd import rmsd
from spyrmsd.adapters.xyzgraph import adjacency_matrix_from_atomic_coordinates


class _FakeGraph:
    def __init__(self, edges):
        self._edges = edges

    def edges(self):
        return list(self._edges)


@pytest.fixture
def fake_xyzgraph(monkeypatch):
    calls = []

    def build_graph(atoms):
        calls.append(
            {
                "atoms": atoms,
            }
        )

        if len(atoms) == 2:
            return _FakeGraph([(0, 1)])

        return _FakeGraph([(0, 1), (0, 2)])

    module = types.SimpleNamespace(build_graph=build_graph)
    monkeypatch.setitem(sys.modules, "xyzgraph", module)
    return calls


def test_xyzgraph_adjacency_returns_symmetric_matrix(fake_xyzgraph):
    atomicnums = np.array([8, 1, 1])
    coordinates = np.array([[0.0, 0.0, 0.0], [0.95, 0.0, 0.0], [-0.24, 0.92, 0.0]])

    adjacency = adjacency_matrix_from_atomic_coordinates(atomicnums, coordinates)

    assert np.array_equal(
        adjacency,
        np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]]),
    )


def test_xyzgraph_adjacency_preserves_input_atom_order(fake_xyzgraph):
    atomicnums = np.array([6, 8, 1])
    coordinates = np.array([[3.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])

    adjacency_matrix_from_atomic_coordinates(atomicnums, coordinates)

    assert fake_xyzgraph[0]["atoms"] == [
        ("C", (3.0, 0.0, 0.0)),
        ("O", (1.0, 0.0, 0.0)),
        ("H", (2.0, 0.0, 0.0)),
    ]


def test_xyzgraph_adjacency_import_error(monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "xyzgraph":
            raise ImportError("No module named 'xyzgraph'")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.delitem(sys.modules, "xyzgraph", raising=False)
    monkeypatch.setattr(builtins, "__import__", fake_import)
    atomicnums = np.array([1, 1])
    coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]])

    with pytest.raises(ImportError, match="xyzgraph"):
        adjacency_matrix_from_atomic_coordinates(atomicnums, coordinates)


def test_xyzgraph_adjacency_rejects_unknown_atomic_number(fake_xyzgraph):
    with pytest.raises(ValueError, match="Unsupported atomic number"):
        adjacency_matrix_from_atomic_coordinates(
            np.array([0]), np.array([[0.0, 0.0, 0.0]])
        )


def test_xyzgraph_adjacency_rejects_out_of_range_edges(monkeypatch):
    module = types.SimpleNamespace(build_graph=lambda atoms: _FakeGraph([(0, 2)]))
    monkeypatch.setitem(sys.modules, "xyzgraph", module)

    with pytest.raises(ValueError, match="out-of-range node indices"):
        adjacency_matrix_from_atomic_coordinates(
            np.array([1, 1]),
            np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]]),
        )


def test_xyzgraph_adjacency_integrates_with_symmrmsd(fake_xyzgraph):
    atomicnums = np.array([8, 1, 1])
    coords_ref = np.array([[0.0, 0.0, 0.0], [0.95, 0.0, 0.0], [-0.24, 0.92, 0.0]])
    coords_pose = coords_ref.copy()

    adjacency = adjacency_matrix_from_atomic_coordinates(atomicnums, coords_ref)

    value = rmsd.symmrmsd(
        coords_ref,
        coords_pose,
        atomicnums,
        atomicnums,
        adjacency,
        adjacency,
    )

    assert value == pytest.approx(0.0, abs=1e-12)


def test_xyzgraph_adjacency_live_smoke():
    try:
        __import__("xyzgraph")
    except Exception:
        pytest.skip("xyzgraph is not importable in this environment")

    atomicnums = np.array([1, 1])
    coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]])

    adjacency = adjacency_matrix_from_atomic_coordinates(atomicnums, coordinates)

    assert adjacency.shape == (2, 2)
    assert np.array_equal(adjacency, adjacency.T)
    assert adjacency[0, 1] == 1

    value = rmsd.symmrmsd(
        coordinates,
        coordinates.copy(),
        atomicnums,
        atomicnums,
        adjacency,
        adjacency,
    )

    assert value == pytest.approx(0.0, abs=1e-12)
