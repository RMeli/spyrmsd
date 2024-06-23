"""
Control skipping of tests according to command line option.

Soource:
    https://docs.pytest.org/en/latest/example/simple.html
"""

import os
import warnings
from collections import namedtuple

import numpy as np
import pytest

import spyrmsd
from spyrmsd import io

Mol = namedtuple("Mol", ["mol", "name", "n_atoms", "n_bonds", "n_h"])


def pytest_addoption(parser):
    parser.addoption(
        "--benchmark", action="store_true", default=False, help="run benchmark"
    )

    parser.addoption(
        "--n-tests",
        default=0,
        type=int,
        help="run n randomly selected tests from large dataset",
    )


def pytest_configure(config):
    """
    Register new markers programmatically:
    - :code:`benchmark`
    - :code:`large`

    Use :code:`pytest -markers` to list all markers.
    """
    config.addinivalue_line("markers", "benchmark: mark test as benchmark")
    config.addinivalue_line("markers", "large: run additional randomly selected tests")

    # Number of system in large test dataset
    pytest.n_systems = 4554


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--benchmark"):
        skip_benchmark = pytest.mark.skip(reason="need --benchmark option to run")

        for item in items:
            if "benchmark" in item.keywords:
                item.add_marker(skip_benchmark)

    if not config.getoption("--n-tests"):
        skip_large = pytest.mark.skip(reason="need --n-tests option to run")

        for item in items:
            if "large" in item.keywords:
                item.add_marker(skip_large)


def pytest_generate_tests(metafunc):
    """
    Parametrise functions based on command line arguments.

    https://docs.pytest.org/en/stable/example/parametrize.html#parametrizing-tests

    Notes
    -----

    The function :code:`test_large/test_rmsd` takes :code:`idx` as argument, which
    need to be parametrised using :code:`--n-tests` from the command line argument.
    """
    if "idx" in metafunc.fixturenames:  # idx is a parameter of test_large/test_rmsd
        n = metafunc.config.getoption("--n-tests")

        metafunc.parametrize("idx", np.random.randint(0, pytest.n_systems, size=n))


@pytest.fixture(autouse=True, params=spyrmsd.available_backends)
def set_backend(request):
    # Capture warning when trying to switch to the same backend
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spyrmsd.set_backend(request.param)


@pytest.fixture(scope="session")
def molpath():
    fdir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(fdir, f"data{os.sep}molecules")


@pytest.fixture
def benzene(molpath):
    mol = io.loadmol(os.path.join(molpath, "benzene.sdf"))
    return Mol(mol, "benzene", 12, 12, 6)


@pytest.fixture
def pyridine(molpath):
    mol = io.loadmol(os.path.join(molpath, "pyridine.sdf"))
    return Mol(mol, "pyridine", 11, 11, 5)


@pytest.fixture
def muparfostat(molpath):
    mol = io.loadmol(os.path.join(molpath, "muparfostat.sdf"))
    return Mol(mol, "muparfostat", 78, 80, 33)


@pytest.fixture(
    params=[
        # (name, n_atoms, n_bonds, n_h)
        ("benzene", 12, 12, 6),
        ("ethanol", 9, 8, 6),
        ("pyridine", 11, 11, 5),
        ("dialanine", 23, 22, 12),
    ]
)
def mol(request, molpath):
    """
    Load molecule as sPyRMSD molecule.
    """

    name, n_atoms, n_bonds, n_h = request.param

    mol = io.loadmol(os.path.join(molpath, f"{name}.sdf"))

    return Mol(mol, name, n_atoms, n_bonds, n_h)


@pytest.fixture
def rawmol(mol, molpath):
    """
    Load molecule as a molecule of the current molecular I/O library.
    """

    RawMol = namedtuple(
        "RawMol", ["mol", "rawmol", "name", "n_atoms", "n_bonds", "n_h"]
    )

    rawmol = io.load(os.path.join(molpath, f"{mol.name}.sdf"))

    return RawMol(mol.mol, rawmol, mol.name, mol.n_atoms, mol.n_bonds, mol.n_h)


@pytest.fixture
def trps(molpath):
    trp_list = []
    for i in range(6):
        trp_list.append(io.loadmol(os.path.join(molpath, f"trp{i}.pdb")))

    return trp_list


@pytest.fixture
def docking_2viz(molpath):
    mols = {}  # Dictionary (pose, molecule)
    for i in [1, 2, 3]:
        mols[i] = io.loadmol(os.path.join(molpath, f"2viz_{i}.sdf"))

    return mols


@pytest.fixture
def docking_1cbr(molpath):
    return [
        io.loadmol(os.path.join(molpath, "1cbr_ligand.mol2")),
        *io.loadallmols(os.path.join(molpath, "1cbr_docking.sdf")),
    ]
