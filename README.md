# sPyRMSD

[![Travis Build Status](https://travis-ci.org/RMeli/spyrmsd.svg?branch=develop)](https://travis-ci.org/RMeli/spyrmsd)
[![Build status](https://ci.appveyor.com/api/projects/status/31ubs980idhv1qw8/branch/develop?svg=true)](https://ci.appveyor.com/project/RMeli/spyrmsd/branch/master)
![flake8](https://github.com/RMeli/spyrmsd/workflows/flake8/badge.svg)
![mypy](https://github.com/RMeli/spyrmsd/workflows/mypy/badge.svg)
[![codecov](https://codecov.io/gh/RMeli/spyrmsd/branch/develop/graph/badge.svg)](https://codecov.io/gh/RMeli/spyrmsd/branch/master)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/RMeli/spyrmsd.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/RMeli/spyrmsd/context:python)

[![Docs](https://img.shields.io/badge/docs-spyrmsd.readthedocs.io-blueviolet)](https://spyrmsd.readthedocs.io)
[![Documentation Status](https://readthedocs.org/projects/spyrmsd/badge/?version=develop)](https://spyrmsd.readthedocs.io/en/develop/?badge=develop)

[![License](https://img.shields.io/github/license/RMeli/pyrmsd?color=%2333BBFF)](https://opensource.org/licenses/MIT)
[![PyPI](https://img.shields.io/badge/PyPI-v0.3.5%20-ff69b4)](https://pypi.org/project/spyrmsd/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/spyrmsd.svg)](https://anaconda.org/conda-forge/spyrmsd)

[![J. Cheminform.](https://img.shields.io/badge/J.%20Cheminform.-10.1186%2Fs13321--020--00455--2-blue)](https://doi.org/10.1186/s13321-020-00455-2)
[![Zenodo](https://zenodo.org/badge/214157073.svg)](https://zenodo.org/badge/latestdoi/214157073)

Python tool for symmetry-corrected RMSD calculations.

---

If you find `spyrmsd` useful, please consider citing the following paper:

```
@article{spyrmsd2020,
  Author = {Meli, Rocco and Biggin, Philip C.},
  Journal = {Journal of Cheminformatics},
  Number = {1},
  Pages = {49},
  Title = {spyrmsd: symmetry-corrected RMSD calculations in Python},
  Volume = {12},
  Year = {2020}
}
```

## Installation

`spyrmsd` is available on [PyPI](https://pypi.org/project/spyrmsd/) and [conda-forge](https://github.com/conda-forge/spyrmsd-feedstock) and can be easily installed from source. See [Dependencies](###Dependencies) for a description of all the dependencies.

_Note_: `spyrmsd` will install [NetworkX](https://networkx.github.io/) (multi-platform). You can install [graph-tool](https://graph-tool.skewed.de/) on macOS and Linux for higher performance.

### PyPI

```bash
pip install spyrmsd
```

### conda

```bash
conda install spyrmsd -c conda-forge
```

### GitHub

```bash
git clone https://github.com/RMeli/spyrmsd.git
cd spyrmsd
pip install .
```

### Dependencies

`spyrmsd` can be used both as a module or as a standalone tool.

#### Module

The following packages are required to use `spyrmsd` as a module:

* [graph-tool](https://graph-tool.skewed.de/) or [NetworkX](https://networkx.github.io/)
* [numpy](https://numpy.org/)
* [scipy](https://www.scipy.org/)
* [qcelemental](http://docs.qcarchive.molssi.org/projects/qcelemental/en/latest/)

_Note_: `spyrmsd` uses [graph-tool](https://graph-tool.skewed.de/) by default but will fall back to [NetworkX](https://networkx.github.io/) if the former is not installed (e.g. on Windows). However, in order to support cross-platform installation [NetworkX](https://networkx.github.io/) is installed by default.

#### Standalone Tool

Additionally, one of the following packages is required to use `spyrmsd` as a standalone tool:

* [Open Babel](http://openbabel.org/)
* [RDKit](https://rdkit.org/)

_Note_: [RDKit](https://rdkit.org/) is not available on PyPI ([Why the RDKit isn't available on PyPi](https://rdkit.blogspot.com/2019/11/why-rdkit-isnt-available-on-pypi.html)). See [RDKit Installation](http://www.rdkit.org/docs/Install.html) for installation instructions.

## Usage

### Standalone

```bash
python -m spyrmsd.spyrmsd -h
```

```text
usage: spyrmsd.py [-h] [-m] [-c] [--hydrogens] [-n] reference molecules [molecules ...]

Python RMSD tool.

positional arguments:
  reference       Reference file
  molecules       Input file(s)

optional arguments:
  -h, --help      show this help message and exit
  -m, --minimize  Minimize (fit)
  -c, --center    Center molecules at origin
  --hydrogens     Keep hydrogen atoms
  -n, --nosymm    No graph isomorphism
```

### Module

```python
from spyrmsd import rmsd
```

#### RMSD

The function  `rmsd.rmsd` computes RMSD without symmetry correction. The atoms are expected to be in the same order for both molecules being compared (no atom matching is performed).

#### Symmetry-Corrected RMSD

The function `rmsd.symmrmsd` computes symmetry-corrected RMSD using molecular graph isomorphisms. Symmetry correction requires molecular adjacency matrices describing the connectivity but needs not the atoms to be in the same order.

Atom matching is performed according to the molecular graph. Therefore, this function should be used when atoms in the molecules being compared are not in the same order (even if there is not symmetry).

## Development

To ensure code quality and consistency the following tools are used during development:

* [black](https://black.readthedocs.io/en/stable/)
* [Flake 8](http://flake8.pycqa.org/en/latest/) (CI)
* [mypy](http://mypy-lang.org/) (CI)

## Copyright

Copyright (c) 2019-2020, Rocco Meli.

## References

References are tracked with [duecredit](https://github.com/duecredit/duecredit/). Run the `credits.sh` script in order to print up-to-date references.

### Acknowledgements

Project based on the [Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version `1.1`.
