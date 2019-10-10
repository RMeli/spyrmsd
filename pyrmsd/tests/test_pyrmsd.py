"""
Unit and regression test for the pyrmsd package.
"""

# Import package, test suite, and other packages as needed
import pyrmsd
import pytest
import sys

def test_pyrmsd_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "pyrmsd" in sys.modules
