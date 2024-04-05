
# `spyrmsd` CHANGELOG

------------------------------------------------------------------------------

## Version 0.7.0

Date:            05/04/2024
Contributors:    @RMeli, @takluyver, @Jnelen

### Added

* Support for `rustworkx` graph library [PR 111 | @RMeli]
* Functionality to manually select the backend from CLI [PR #108 | @RMeli]
* Functionality to manually select the backend [PR  #107 | @Jnelen]
* Python `3.12` to CI [PR  #102 | @RMeli]
* macOS M1 (`macoOS-14`) to CI [PR  #102 | @RMeli]

### Changed

* Molecular graphs cache to cache by backend [PR  #107 | @Jnelen]
* Python build system to use flit_core directly [PR #103 | @takluyver]
* Minimum version of Python to `3.9` (to reduce CI matrix) [PR  #102 | @RMeli]

### Fixed

* Failing tests with `pytest=8.0.0` [PR #101 | @RMeli]

### Improved

* Messages for `NotImplementedError` exceptions [PR #90 | @RMeli]

### Removed

* `.gitattributes` and `.lgtm.yaml` stale files
* `versioneer` [PR #91 | @RMeli]

## Version 0.6.0

Date:            08/09/2023
Contributors:    @RMeli

### Improved

* Documentation on loading RDKit and Open Babel molecules [PR #77 | @RMeli]

### Changed

* Versions of `pre-commit` hooks tools [PR #85 | @RMeli]
* Version of several GitHub actions [PR #76 | @RMeli]
* Versioneer to `0.28` [PR #76 | @RMeli]
* Minimum Python version to `3.7` to match CI [PR #76 | @RMeli]
* Code according to `black==23.3.0` [PR #76 | @RMeli]

### Added

* `pre-commit` CI action [PR #85 | @RMeli]
* `extras_require` to `setup.py` for RDKit and Open Babel [PR #84 | @RMeli]
* Error message when `spyrmsd` is used as module but neither OpenBabel nor RDKit are installed [PR #81 | @RMeli]
* Pin to `sphinx<7` to avoid compatibility with RTD theme [PR #77 | @RMeli]

### Removed

* `ubuntu-latest-3.7-rdkit-gt` configuration from CI due to continuous failures [PR #84 | @RMeli]
* Outdated information about RDKit from the documentation [PR #84 | @RMeli]
* Support for Open Babel 2 [PR #84 | @RMeli]
* LGTM badge and code annotations [PR #76 | @RMeli]

------------------------------------------------------------------------------

## Version 0.5.2

Date:            23/02/2022
Contributors:    @RMeli

### Fixed

* Inconsistent number of nodes with `graph-tool` for disconnected graphs [PR #61 | @RMeli]

### Improved

* Support for more types of node properties (including strings) with `graph-tool` [PR #64 | @RMeli]

### Changed

* `ValueError` exception into `NonIsomorphicGraphs(ValueError)` exception [PR #65 | @RMeli]

### Added

* Warning for disconnected graphs [PR #61| @RMeli]

------------------------------------------------------------------------------

## Version 0.5.1

Date:            21/09/2021
Contributors:    @RMeli

### Fixed

* Wrong covalent radius in `graph.adjacency_matrix_from_atomic_coordinates()` [PR #58 | @RMeli]

### Added

* [pre-commit](https://pre-commit.com/) configuration file [PR #57 | @RMeli]
* Support for gzip-compressed files (`.gz`) [PR #56 | @RMeli]

### Removed

* Dependency `QCElemental` [PR #58 | @RMeli]

------------------------------------------------------------------------------

## Version 0.5.0

Date:            21/06/2021

Contributors:    @RMeli

### Added

* `molecule.Molecule` constructor from RDKit molecule [PR #50 | @RMeli]
* `molecule.Molecule` constructor from Open Babel molecule [PR #50 | @RMeli]
* `--n-tests` option for `pytest` [PR #44 | @RMeli]

### Improved

* `spyrmsd.rmsdwrapper` to deal with single molecule [PR #51 | @RMeli]
* Issue template [PR #46 | @RMeli]
* Speed of computation of squared pairwise distances [PR #45 | @RMeli]

### Changed

* `spyrmsd` standalone tool now invoked with `python -m spyrmsd -h` [PR #52 | @RMeli]
* Moved `spyrmsd.coords_from_molecule` to `molecule` module [PR #52 | @RMeli]
* Moved `spyrmsd.rmsdwrapper` to `rmsd` module [PR #52 | @RMeli]
* Long tests no longer run in CI [PR #44 | @RMeli]

### Removed

* `spyrmsd` module [PR #52 | @RMeli]
* Travis CI and AppVeyor bindings [PR #44 | @RMeli]
* `--long` option for `pytest` [PR #44 | @RMeli]

------------------------------------------------------------------------------
