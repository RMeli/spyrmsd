"""
spyrmsd

Python tool for symmetry-corrected RMSD
"""

from setuptools import find_packages, setup

import versioneer

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
    author_email="rocco.meli@biodtp.ox.ac.uk",
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    url="https://spyrmsd.readthedocs.io",
    install_requires=["numpy", "scipy", "qcelemental", "networkx>=2"],
    extras_require={"bib": ["duecredit"]},
    platforms=["Linux", "Mac OS-X", "Unix", "Windows"],
    python_requires=">=3.6",
)
