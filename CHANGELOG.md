
# `spyrmsd` CHANGELOG

------------------------------------------------------------------------------

## Version X.Y.Z

Date:            DD/MM/YYYY
Contributors:    @RMeli

### Added

* Added support for gzip-compressed files (`.gz`) [PR #56 | @RMeli]

## Version 0.5.0

Date:            21/06/2021

Contributors:    @RMeli

### Added

* Added `molecule.Molecule` constructor from RDKit molecule [PR #50 | @RMeli]
* Added `molecule.Molecule` constructor from Open Babel molecule [PR #50 | @RMeli]
* Added `--n-tests` option for `pytest` [PR #44 | @RMeli]

### Improved

* Improved `spyrmsd.rmsdwrapper` to deal with single molecule [PR #51 | @RMeli]
* Improved issue template [PR #46 | @RMeli]
* Improved speed of computation of squared pairwise distances [PR #45 | @RMeli]

### Changed

* `spyrmsd` standalone tool now invoked with `python -m spyrmsd -h` [PR #52 | @RMeli]
* Moved `spyrmsd.coords_from_molecule` to `molecule` module [PR #52 | @RMeli]
* Moved `spyrmsd.rmsdwrapper` to `rmsd` module [PR #52 | @RMeli]
* Long tests no longer run in CI [PR #44 | @RMeli]

### Removed

* Removed `spyrmsd` module [PR #52 | @RMeli]
* Removed Travis CI and AppVeyor bindings [PR #44 | @RMeli]
* Removed `--long` option for `pytest` [PR #44 | @RMeli]

------------------------------------------------------------------------------
