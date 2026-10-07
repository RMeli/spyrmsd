import builtins
import sys
import types

import numpy as np
import pytest

from spyrmsd import __main__ as cli
from spyrmsd import graph
from spyrmsd.adapters import xyzgraph as xyzgraph_adapter


@pytest.fixture(autouse=True)
def reset_graph_builder():
    graph.set_graph_builder("simple")
    yield
    graph.set_graph_builder("simple")


def run_cli(monkeypatch, *args):
    monkeypatch.setattr(sys, "argv", ["spyrmsd", *args])
    return cli.main()


def test_main_requires_rdkit(monkeypatch):
    monkeypatch.setattr(cli.importlib.util, "find_spec", lambda name: None)

    with pytest.raises(ImportError, match="RDKit not found"):
        run_cli(monkeypatch, "ref.sdf", "poses.sdf")


def test_main_reports_missing_reference(monkeypatch, capsys):
    monkeypatch.setattr(cli.importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(
        cli.io, "loadmol", lambda path: (_ for _ in ()).throw(OSError())
    )

    assert run_cli(monkeypatch, "missing.sdf", "poses.sdf") == -1
    assert capsys.readouterr().err == "ERROR: Reference file not found.\n"


def test_main_reports_missing_molecules(monkeypatch, capsys):
    fake_mol = types.SimpleNamespace(
        atomicnums=np.array([1]),
        coordinates=np.zeros((1, 3)),
    )
    monkeypatch.setattr(cli.importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(cli.io, "loadmol", lambda path: fake_mol)
    monkeypatch.setattr(
        cli.io, "loadallmols", lambda path: (_ for _ in ()).throw(OSError())
    )

    assert run_cli(monkeypatch, "ref.sdf", "missing.sdf") == -1
    assert capsys.readouterr().err == "ERROR: Molecule file(s) not found.\n"


def test_main_default_path_uses_io(monkeypatch, capsys):
    fake_mol = types.SimpleNamespace(
        atomicnums=np.array([8, 1, 1]),
        coordinates=np.zeros((3, 3)),
        adjacency_matrix=np.zeros((3, 3), dtype=int),
    )

    monkeypatch.setattr(
        cli.importlib.util,
        "find_spec",
        lambda name: object() if name == "rdkit" else None,
    )
    monkeypatch.setattr(cli.io, "loadmol", lambda p: fake_mol)
    monkeypatch.setattr(cli.io, "loadallmols", lambda p: [fake_mol, fake_mol])
    monkeypatch.setattr(
        cli,
        "rmsdwrapper",
        lambda ref, mols, symmetry, center, minimize, strip, cache=True: [0.1, 0.2],
    )

    rc = run_cli(monkeypatch, "ref.sdf", "poses.sdf")

    assert rc == 0
    assert graph.get_graph_builder() == "simple"
    assert capsys.readouterr().out.strip().splitlines() == ["0.10000", "0.20000"]


def test_main_sets_graph_backend(monkeypatch):
    fake_mol = types.SimpleNamespace(
        atomicnums=np.array([1]),
        coordinates=np.zeros((1, 3)),
    )
    selected = []

    monkeypatch.setattr(cli.importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(cli.io, "loadmol", lambda path: fake_mol)
    monkeypatch.setattr(cli.io, "loadallmols", lambda path: [fake_mol])
    monkeypatch.setattr(cli.spyrmsd, "set_backend", selected.append)
    monkeypatch.setattr(cli, "rmsdwrapper", lambda *args, **kwargs: [0.0])

    assert (
        run_cli(monkeypatch, "--graph-backend", "networkx", "ref.sdf", "poses.sdf") == 0
    )
    assert selected == ["networkx"]


def test_main_xyzgraph_builder_rebuilds_adjacency(monkeypatch, capsys):
    fake_mol = types.SimpleNamespace(
        atomicnums=np.array([8, 1, 1]),
        coordinates=np.zeros((3, 3)),
        adjacency_matrix=np.zeros((3, 3), dtype=int),
    )
    monkeypatch.setattr(
        xyzgraph_adapter,
        "adjacency_matrix_from_atomic_coordinates",
        lambda aprops, coordinates: np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]),
    )

    monkeypatch.setattr(
        cli.importlib.util,
        "find_spec",
        lambda name: object() if name == "rdkit" else None,
    )
    monkeypatch.setattr(cli.io, "loadmol", lambda p: fake_mol)
    monkeypatch.setattr(cli.io, "loadallmols", lambda p: [fake_mol])
    seen = {}

    def fake_rmsdwrapper(ref, mols, symmetry, center, minimize, strip, cache=True):
        seen["ref"] = ref
        seen["mols"] = mols
        return [0.0]

    monkeypatch.setattr(cli, "rmsdwrapper", fake_rmsdwrapper)

    rc = run_cli(monkeypatch, "--graph-builder", "xyzgraph", "ref.xyz", "poses.xyz")

    assert rc == 0
    assert graph.get_graph_builder() == "xyzgraph"
    assert np.array_equal(
        seen["ref"].adjacency_matrix,
        np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]),
    )
    assert np.array_equal(
        seen["mols"][0].adjacency_matrix,
        np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]),
    )
    assert capsys.readouterr().out.strip() == "0.00000"


def test_main_xyzgraph_builder_missing_dependency(monkeypatch, capsys):
    fake_mol = types.SimpleNamespace(
        atomicnums=np.array([8, 1, 1]),
        coordinates=np.zeros((3, 3)),
        adjacency_matrix=np.zeros((3, 3), dtype=int),
    )
    real_import = builtins.__import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "xyzgraph":
            raise ImportError("No module named 'xyzgraph'")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(
        cli.importlib.util,
        "find_spec",
        lambda name: object() if name == "rdkit" else None,
    )
    monkeypatch.setattr(cli.io, "loadmol", lambda p: fake_mol)
    monkeypatch.setattr(cli.io, "loadallmols", lambda p: [fake_mol])
    monkeypatch.setattr(builtins, "__import__", fake_import)

    rc = run_cli(monkeypatch, "--graph-builder", "xyzgraph", "ref.xyz", "poses.xyz")

    assert rc == -1
    err = capsys.readouterr().err
    assert err.startswith("ERROR: xyzgraph is required for the xyzgraph graph builder")


def test_main_verbose_prints_library_and_builder(monkeypatch, capsys):
    fake_mol = types.SimpleNamespace(
        atomicnums=np.array([8, 1, 1]),
        coordinates=np.zeros((3, 3)),
        adjacency_matrix=np.zeros((3, 3), dtype=int),
    )

    monkeypatch.setattr(
        cli.importlib.util,
        "find_spec",
        lambda name: object() if name == "rdkit" else None,
    )
    monkeypatch.setattr(cli.io, "loadmol", lambda p: fake_mol)
    monkeypatch.setattr(cli.io, "loadallmols", lambda p: [fake_mol])
    monkeypatch.setattr(
        cli,
        "rmsdwrapper",
        lambda ref, mols, symmetry, center, minimize, strip, cache=True: [0.0],
    )

    rc = run_cli(monkeypatch, "--verbose", "ref.sdf", "poses.sdf")

    assert rc == 0
    lines = capsys.readouterr().out.strip().splitlines()
    assert lines[0].startswith("Graph library:")
    assert lines[1] == "Graph builder: simple"
    assert lines[2] == "0.00000"
