
# `spyrmsd` CHANGELOG

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
