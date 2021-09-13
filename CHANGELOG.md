# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Support for xtb optimizations, single-point energies, and molecular dynamics. (xtb)
- Support for additional functionals. (ORCA)
- Support for coupled cluster calculations. (ORCA)
- Track Hartree&ndash;Fock type for ab initio methods. (ORCA)
- A script to combine JSON files recursively.
- Adds `n_electrons` to JSON files. (ORCA)

### Changed

- Fixed handling of errors in recursive QCJSON creator operations.
- Fixed including and excluding keywords in paths.
- More consistent information about frozen core electrons. Also changed
  `'frozencore'` to ``'frozen_core'``.
- Parsing errors will not stop qcjson-creator.py if debug is False.
- UHF spin contamination correctly parsed for DFT methods. (ORCA)
- Adding missing &omega;B97M-V and &omega;B97M-D3BJ functionals. (ORCA)

## [0.2.0] 2021-02-09

### Added

- Frequency, normal modes, and thermochemistry data. (ORCA)
- Mulliken and Loewdin charges extraction. (ORCA)
- Dipole moment property. (ORCA)

### Changed

- Properly keeps track of auxillary basis sets. (ORCA)
- Improved exclusion and inclusions system for filtering files.

## [0.1.1] - 2021-01-06

### Fixed

- cclib version requirement.
- Recursive option would save in current directory and not in the same directory of the output file.
- get_json would incorrectly catch KeyboardInterrupt exception.

## [0.1.0] - 2021-01-05

### Added

- Custom parser for ORCA information such as integration grid, scf energy contributions (e.g., one-electron and two-electron energies), MP2 correlation energies, RI approximations. (ORCA)
- Debug option to raise errors instead of skipping over files.
- Alpha and beta electron HOMO and LUMO information. (ORCA)
- 'return_energy' property regardless of driver.

### Changed

- Nest iterations into a list instead of having int labels.
- Standardized getting SCF, MP, and CC energies from cclib.
- Requires outfile path to initialize json classes.
- Write each JSON file directly after parsing instead of all at the end. That
  way if the script crashes the proceeding JSON files are already written.

## [0.0.1] - 2021-01-03

- Initial release!
