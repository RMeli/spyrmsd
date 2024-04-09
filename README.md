# sPyRMSD

[![pytest](https://github.com/RMeli/spyrmsd/actions/workflows/pytest.yml/badge.svg?branch=develop)](https://github.com/RMeli/spyrmsd/actions/workflows/pytest.yml)
![flake8](https://github.com/RMeli/spyrmsd/workflows/flake8/badge.svg)
![mypy](https://github.com/RMeli/spyrmsd/workflows/mypy/badge.svg)
[![codecov](https://codecov.io/gh/RMeli/spyrmsd/branch/develop/graph/badge.svg)](https://codecov.io/gh/RMeli/spyrmsd/branch/master)

[![Docs](https://img.shields.io/badge/docs-spyrmsd.readthedocs.io-blueviolet)](https://spyrmsd.readthedocs.io)
[![Documentation Status](https://readthedocs.org/projects/spyrmsd/badge/?version=develop)](https://spyrmsd.readthedocs.io/en/develop/?badge=develop)

[![License](https://img.shields.io/github/license/RMeli/pyrmsd?color=%2333BBFF)](https://opensource.org/licenses/MIT)
[![PyPI](https://img.shields.io/badge/PyPI-v0.7.0%20-ff69b4)](https://pypi.org/project/spyrmsd/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/spyrmsd.svg)](https://anaconda.org/conda-forge/spyrmsd)

[![J. Cheminform.](https://img.shields.io/badge/J.%20Cheminform.-10.1186%2Fs13321--020--00455--2-blue)](https://doi.org/10.1186/s13321-020-00455-2)
[![Zenodo](https://zenodo.org/badge/214157073.svg)](https://zenodo.org/badge/latestdoi/214157073)

Python tool for symmetry-corrected RMSD calculations.

---

If you find `spyrmsd` useful, please consider citing the following paper:

```text
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

> [!NOTE]
> `spyrmsd` will install [NetworkX](https://networkx.github.io/) (multi-platform). You can install the other backends for higher performance.

> [!WARNING]
> If `spyrmsd` is used as a standalone tool, it is required to install either [RDKit](https://rdkit.org/) or [Open Babel](http://openbabel.org/). Neither is automatically installed with `pip` nor `conda`.

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

* [numpy](https://numpy.org/)
* [scipy](https://www.scipy.org/)

One of the following graph libraries is required:
* [graph-tool]
* [NetworkX]
* [rustworkx]

> [!NOTE]
> `spyrmsd` uses the following priority when multiple graph libraries are present: [graph-tool], [NetworkX], [rustworkx]. *This order might change. Use `set_backend` to ensure you are always using the same backend, if needed.* However, in order to support cross-platform installation [NetworkX](https://networkx.github.io/) is installed by default, and the other graph library need to be installed manually.

#### Standalone Tool

Additionally, one of the following packages is required to use `spyrmsd` as a standalone tool:

* [Open Babel](http://openbabel.org/)
* [RDKit](https://rdkit.org/)

## Usage

### Standalone Tool

`spyrmsd` provides a convenient CLI tool. See `spyrmsd`'s `--help` for the usage:

```bash
python -m spyrmsd -h
```

### Module

```python
from spyrmsd import rmsd
```

#### RMSD

The function  `rmsd.rmsd` computes RMSD without symmetry correction. The atoms are expected to be in the same order for both molecules being compared (no atom matching is performed).

```python
def rmsd(
    coords1: np.ndarray,    # Coordinates of molecule 1
    coords2: np.ndarray,    # Coordinates of molecule 2
    aprops1: np.ndarray,    # Atomic properties of molecule 1
    aprops2: np.ndarray,    # Atomic properties of molecule 2
    center: bool = False,   # Flag to center molecules at origin
    minimize: bool = False, # Flag to compute minimum RMSD
    atol: float = 1e-9,     # Numerical tolerance for QCP method
)
```

> [!NOTE]
> Atomic properties (`aprops`) can be any Python object when using [NetworkX] and [rustworkx], or integers, floats, or strings when using [graph-tool].

#### Symmetry-Corrected RMSD

The function `rmsd.symmrmsd` computes symmetry-corrected RMSD using molecular graph isomorphisms. Symmetry correction requires molecular adjacency matrices describing the connectivity but needs not the atoms to be in the same order.

Atom matching is performed according to the molecular graph. This function should also be used when atoms in the molecules being compared are not in the same order (even if there is not symmetry to be accounted for).

```python
def symmrmsd(
    coordsref: np.ndarray,                       # Reference coordinated
    coords: Union[np.ndarray, List[np.ndarray]], # Coordinates (one set or multiple sets)
    apropsref: np.ndarray,                       # Reference atomic properties
    aprops: np.ndarray,                          # Atomic properties
    amref: np.ndarray,                           # Reference adjacency matrix
    am: np.ndarray,                              # Adjacency matrix
    center: bool = False,                        # Flag to center molecules at origin
    minimize: bool = False,                      # Flag to compute minimum RMSD
    cache: bool = True,                          # Cache graph isomorphisms
    atol: float = 1e-9,                          # Numerical tolerance for QCP method
)
```

> [!NOTE]
> Atomic properties (`aprops`) can be any Python object when using [NetworkX] and [rustworkx], or integers, floats, or strings when using [graph-tool](https://graph-tool.skewed.de/).

#### Select Backend

`spyrmsd` supports the following graph libraries for the calculation of graph isomorphisms:
* [graph-tool]
* [NetworkX]
* [rustworkx]

 You can check which backend is being used with

```python
spyrmsd.get_backend()
```

You can also manually select your preferred backend with

```python
spyrmsd.set_backend("networkx")
# spyrmsd uses NetworkX
spyrmsd.set_backend("graph_tool")
# spyrmsd uses graph_tool
```

The available backends (which depend on the installed dependencies) are stored in `spyrmsd.available_backends`.

## Development

To ensure code quality and consistency the following tools are used during development:

* [black](https://black.readthedocs.io/en/stable/)
* [Flake 8](http://flake8.pycqa.org/en/latest/) (CI)
* [isort]()
* [mypy](http://mypy-lang.org/) (CI)

Pre-commit `git` hooks can be installed with [pre-commit](https://pre-commit.com/).

## Copyright

Copyright (c) 2019-2024, Rocco Meli.

## References

References are tracked with [duecredit](https://github.com/duecredit/duecredit/). Run the `credits.sh` script in order to print up-to-date references.

### Acknowledgements

Project based on the [Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version `1.1`.

[rustworkx]: https://www.rustworkx.org
[NetworkX]: https://networkx.github.io/
[graph-tool]: https://graph-tool.skewed.de/
