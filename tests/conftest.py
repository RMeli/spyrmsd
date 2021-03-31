"""
Control skipping of tests according to command line option.

Soource:
    https://docs.pytest.org/en/latest/example/simple.html
"""
import numpy as np
import pytest


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
