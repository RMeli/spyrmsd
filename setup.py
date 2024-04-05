"""
spyrmsd

Python tool for symmetry-corrected RMSD
"""

from setuptools import find_packages, setup

short_description = __doc__.split("\n")

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except Exception:
    long_description = "\n".join(short_description[2:])


setup(
    # Self-descriptive entries which should always be present
    name="spyrmsd",
    author="Rocco Meli",
    author_email="rocco.meli@cscs.ch",
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version="0.7.0-dev",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    url="https://spyrmsd.readthedocs.io",
    install_requires=["numpy", "scipy", "networkx>=2"],
    extras_require={
        "bib": ["duecredit"],
        "rdkit": ["rdkit"],
        "openbabel": ["openbabel"],
    },
    platforms=["Linux", "Mac OS-X", "Unix", "Windows"],
    python_requires=">=3.9",
)
