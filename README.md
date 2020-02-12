# sPyRMSD

[![Travis Build Status](https://travis-ci.org/RMeli/spyrmsd.svg?branch=develop)](https://travis-ci.org/RMeli/spyrmsd)
[![Build status](https://ci.appveyor.com/api/projects/status/31ubs980idhv1qw8/branch/develop?svg=true)](https://ci.appveyor.com/project/RMeli/spyrmsd/branch/master)
[![codecov](https://codecov.io/gh/RMeli/spyrmsd/branch/develop/graph/badge.svg)](https://codecov.io/gh/RMeli/spyrmsd/branch/master)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/RMeli/spyrmsd.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/RMeli/spyrmsd/context:python)
[![Documentation Status](https://readthedocs.org/projects/spyrmsd/badge/?version=develop)](https://spyrmsd.readthedocs.io/en/develop/?badge=develop)

[![DOI](https://zenodo.org/badge/214157073.svg)](https://zenodo.org/badge/latestdoi/214157073)
[![Docs](https://img.shields.io/badge/docs-pyrmsd.readthedocs.io-blueviolet)](https://pyrmsd.readthedocs.io)
[![PyPI](https://img.shields.io/badge/PyPI-0.2.1%20-ff69b4)](https://pypi.org/project/spyrmsd/)
[![License](https://img.shields.io/github/license/RMeli/pyrmsd?color=%2333BBFF)](https://opensource.org/licenses/MIT)

Python-first tool for symmetry-corrected RMSD calculations.

## Installation

### PyPI

```bash
pip install spyrmsd
```

This install `spyrmsd` as a library. In order to use the standalone RMSD calculation tool, you will need to install [Open Babel](http://openbabel.org/) or [RDKit](https://rdkit.org/).

### GitHub

```bash
git clone https://github.com/RMeli/spyrmsd.git
cd spyrmsd
pip install .
```

### conda

## Usage

### Standalone

```bash
python -m spyrmsd.spyrmsd -h
```

```text
usage: spyrmsd.py [-h] [-m] [-c] [--hydrogens] [-n]
                  reference molecules [molecules ...]

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

The function  `rmsd.rmsd` computes RMSD without symmetry correction. The atoms are expected to be in the same order for both molecule being compared (no atom matching is performed).

#### Symmetric RMSD

The function `rmsd.rmsd_isomorphic` computes symmetry-corrected RMSD, using molecular graph isomorphisms. Symmetry correction requires molecular adjacency matrices but needs not the atoms to be in the same order.

## Contributions

### List of Contributors

| Name               | GitHub          |
| :----------------: | :-------------: |
| Rocco Meli         | @RMeli          | 

### Tools

#### Formatting

The code is automatically formatted using [black](https://black.readthedocs.io/en/stable/):

```bash
black .
```

#### Style

Code style is enforced using [Flake 8](http://flake8.pycqa.org/en/latest/)

```bash
flake8
```

#### Static Checks

Static checks are performed using [mypy](http://mypy-lang.org/)

```bash
mypy
```

## Deployment

### PyPI

Build wheel and sdist:

```bash
flit build
```

Upload wheel and sdist on [PyPI](https://pypi.org/):

```bash
flit publish
```

## References

| Method    | Reference                                          | DOI |
| :-------: | -------------------------------------------------- | :--: |
| QCP       | D. L. Theobald, Acta Crys. A**61**, 478-480 (2005) | [![doi](https://img.shields.io/badge/doi-10.1107%2FS0108767305015266-blue)](https://doi.org/10.1107/S0108767305015266) |
| Hungarian | W. J. Allen and R. C. Rizzo, J. Chem. Inf. Model. **54**, 518-529 (2014) | [![doi](https://img.shields.io/badge/doi-10.1021%2Fci400534h-blue)](https://doi.org/10.1021/ci400534h) |
| Connectivity | E. C. Meng and R. A. Lewis, J. Comp. Chem. **12**, 891-898 (1991) | [![doi](https://img.shields.io/badge/doi-10.1002%2Fjcc.540120716-blue)](https://doi.org/10.1002/jcc.540120716) |

## Copyright

Copyright (c) 2019-2020, Rocco Meli.

### Acknowledgements

Project based on the [Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version `1.1`.
