"""
Python RMSD tool with symmetry correction.
"""

from .due import due, Doi
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

# This will print latest Zenodo version
due.cite(
    Doi("10.5281/zenodo.3631876"),  # lgtm[py/procedure-return-value-used]
    path="spyrmsd",
    description="spyRMSD: Symmetry-Corrected RMSD in Python",
)
