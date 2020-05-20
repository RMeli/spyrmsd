"""
Test imported package and modules.

This test is useful to generate a citation report:

.. code

   python -m duecite test_import.py

"""

import sys

import pytest

from spyrmsd import graph, hungarian, io, molecule, qcp, rmsd, spyrmsd, utils


@pytest.mark.parametrize(
    "name",
    [
        # Package
        spyrmsd.__name__,
        # Modules
        graph.__name__,
        hungarian.__name__,
        io.__name__,
        molecule.__name__,
        qcp.__name__,
        rmsd.__name__,
        spyrmsd.__name__,
        utils.__name__,
    ],
)
def test_imported(name):
    assert name in sys.modules
