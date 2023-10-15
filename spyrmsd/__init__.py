"""
Python RMSD tool with symmetry correction.
"""

from .due import Doi, due

__version__ = "0.7.0-dev"

# This will print latest Zenodo version
due.cite(
    Doi("10.5281/zenodo.3631876"),
    path="spyrmsd",
    description="spyrmsd",
)

due.cite(
    Doi("10.1186/s13321-020-00455-2"),
    path="spyrmsd",
    description="spyrmsd",
)
