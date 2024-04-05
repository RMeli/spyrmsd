"""
Python RMSD tool with symmetry correction.
"""

from .due import Doi, due

# Make the backend related functions available from base spyrmsd
# Add noqa to avoid flake8 error
from .graph import _available_backends as available_backends  # noqa: F401
from .graph import _get_backend as get_backend  # noqa: F401
from .graph import _set_backend as set_backend  # noqa: F401

__version__ = "0.7.0"

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
