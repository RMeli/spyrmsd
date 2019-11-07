PyRMSD
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/RMeli/pyrmsd.svg?branch=master)](https://travis-ci.org/RMeli/pyrmsd)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/rhd5wi1ce7i24hgb/branch/master?svg=true)](https://ci.appveyor.com/project/rhd5wi1ce7i24hgb/pyrmsd/branch/master)
[![codecov](https://codecov.io/gh/RMeli/pyrmsd/branch/master/graph/badge.svg)](https://codecov.io/gh/RMeli/pyrmsd/branch/master)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/RMeli/pyrmsd.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/RMeli/pyrmsd/context:python)
[![Documentation Status](https://readthedocs.org/projects/pyrmsd/badge/?version=latest)](https://pyrmsd.readthedocs.io/en/latest/?badge=latest)

[![Docs](https://img.shields.io/badge/docs-pyrmsd.readthedocs.io-blueviolet)](https://pyrmsd.readthedocs.io)
[![License](https://img.shields.io/github/license/RMeli/pyrmsd?color=%2333BBFF)](https://opensource.org/licenses/MIT)

Python-first tool for symmetry-corrected  RMSD calculations.

## Installation

### GitHub

```bash
git clone https://github.com/RMeli/pyrmsd.git
cd pyrmsd
pip install .
```

### pip

### conda

## Usage

### Standalone

```bash
python -m pyrmsd.pyrmsd -h
```

```text
usage: pyrmsd.py [-h] [-s] [-m] [-v] reference molecules [molecules ...]

Python RMSD tool.

positional arguments:
  reference       Reference file
  molecules       Input file

optional arguments:
  -h, --help      show this help message and exit
  -s, --strip     Strip H atoms
  -m, --minimize  Minimize (fit)
  -v, --verbose   Verbose output
```

## Contributions

### Formatting

The code is automatically formatted using [black](https://black.readthedocs.io/en/stable/):

```bash
black .
```

### Style

Code style is enforced using [Flake 8](http://flake8.pycqa.org/en/latest/)

```bash
flake8
```

### Static Checks

Static checks are performed using [mypy](http://mypy-lang.org/)

```bash
mypy
```

## References

| Method    | Reference                                          | DOI |
| :-------: | -------------------------------------------------- | :--: |
| QCP       | D. L. Theobald, Acta Crys. A**61**, 478-480 (2005) | [![doi](https://img.shields.io/badge/doi-10.1107%2FS0108767305015266-blue)](https://doi.org/10.1107/S0108767305015266) |
| Hungarian | W. J. Allen and R. C. Rizzo, J. Chem. Inf. Model. **54**, 518-529 (2014) | [![doi](https://img.shields.io/badge/doi-10.1021%2Fci400534h-blue)](https://doi.org/10.1021/ci400534h)

## Copyright

Copyright (c) 2019, Rocco Meli.

### Acknowledgements
 
Project based on the [Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version `1.1`.
