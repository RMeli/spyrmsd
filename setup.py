"""
spyrmsd

Python tool for symmetry-corrected RMSD
"""
import sys

from setuptools import find_packages, setup

import versioneer

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest-runner"] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except Exception:
    long_description = "\n".join(short_description[2:])


setup(
    # Self-descriptive entries which should always be present
    name="spyrmsd",
    author="Rocco Meli",
    author_email="rocco.meli@biodtp.ox.ac.uk",
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="MIT",
    packages=find_packages(),
    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,
    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,
    url="https://spyrmsd.readthedocs.io",
    # install_requires=[], # Required packages, pulls from pip if needed; do not use for Conda deployment
    platforms=["Linux", "Mac OS-X", "Unix", "Windows"],
    python_requires=">=3.6",
    # zip_safe=False,
)
