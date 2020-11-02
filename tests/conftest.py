"""
Control skipping of tests according to command line option.

Soource:
    https://docs.pytest.org/en/latest/example/simple.html
"""
import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--benchmark", action="store_true", default=False, help="run benchmark"
    )

    parser.addoption(
        "--large",
        action="store_true",
        default=False,
        help="run large number of randomly selected tests",
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "benchmark: mark test as benchmark")
    config.addinivalue_line("markers", "large: mark test as large")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--benchmark"):
        skip_benchmark = pytest.mark.skip(reason="need --benchmark option to run")

        for item in items:
            if "benchmark" in item.keywords:
                item.add_marker(skip_benchmark)

    if not config.getoption("--large"):
        skip_large = pytest.mark.skip(reason="need --large option to run")

        for item in items:
            if "large" in item.keywords:
                item.add_marker(skip_large)
