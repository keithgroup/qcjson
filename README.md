# qcjson

Reproducibility is a cornerstone of science, and computational chemistry is no different.
With decreasing levels of transparency and reproducibility in computational chemistry, researchers should provide source codes, raw output files, input files, or some sort of machine-readable database (e.g., XML and JSON).
Since commercial software is commonly used, this incentivizes providing output files and a machine readable database instead.
This repository contains scripts for the creation and use of Quantum Chemistry JSONs, QCJSONs, for computational chemistry output files using a modified [QCSchema](https://github.com/MolSSI/QCSchema).

## Table of Contents

- [Table of Contents](#table-of-contents)
- [Supported packages](#supported-packages)
- [Requirements](#requirements)
- [Installation](#installation)
- [Scripts](#scripts)
  - [qcjson-creator](#qcjson-creator)
    - [qcjson-creator options](#qcjson-creator-options)
    - [qcjson-creator examples](#qcjson-creator-examples)
- [QCJSON Format](#qcjson-format)
  - [Keys](#keys)
- [License](#license)

## Supported packages

The following computational packages and calculation types are currently supported (versions tested):

- ORCA (v4.2) optimizations, single point energies, gradients, frequencies.

## Requirements

Python 3 with the following packages:

- cclib (>=1.7).

## Installation

1. Clone the repo

   ```text
   git clone https://github.com/keithgroup/qcjson
   ```

2. Install the qcjson package using pip.

   ```text
   cd qcjson
   pip3 install --user .
   ```

3. Add scripts to your path to use them from anywhere (optional).
   For example, if you cloned the repository in your home directory on Linux you would add `export PATH="$HOME/qcjson/scripts:$PATH"` to your ~/.bashrc file.

## Scripts

Python scripts are provided for common uses.

### qcjson-creator

Python script that creates QCJSON files from computational chemistry output files.

#### qcjson-creator options

- `-o`, `--overwrite`: Overwrite JSON file if it already exists.
- `-r`, `--recursive`: Recursively find output files to make QCJSONs.
- `-p`, `--prettify`: Indent JSON properties.
- `--save_dir`: Path to directory to store JSON files; for example, `--save_dir json_dir`.
  Defaults to current directory.
- `--exclude`: Does not include file paths that have these strings in their path; for example, `--exclude bad bp86` will not include any files with 'bad' or 'bp86' in the path.
- `--include`: Will only include files that have all of these strings in their path; for example, `--include b3lyp` will only include files with `b3lyp` in the path.
- `-d`, `--debug`: Will raise errors instead of skipping failed JSONs.

#### qcjson-creator examples

There are multiple example calculations and QCJSONs in the `data` directory.
For example, this is how to create QCJSONs for all of the examples in the `orca` directory.

```text
$ qcjson-creator.py -rop .
QCJSON creator v0.2.0
Written by Alex M. Maldonado (@aalexmmaldonado)
Energies and distances are Hartrees and Angstroms

Looking for output files in ./, recursively
Found 8 output files

Making QCJSON for 4h2o.abc0.iter2-orca.opt-mp2.def2tzvp
Making QCJSON for 12h2o.su.etal-orca.sp-mp2.def2tzvp.frozencore
Making QCJSON for alo3h3.h2o-orca.freq-bp86.def2svp.cpcm
Making QCJSON for neb-bare.0-orca.sp-bp86.def2tzvp.smd
Making QCJSON for 5h2o.abc0.iter1.mp2.md.300k.iter1.mol0,1,2,3,4-orca.engrad-mp2.def2tzvp
Making QCJSON for alo4h5-orca.freq-bp86.def2svp
Making QCJSON for neb-bare.0-orca.sp-wb97x.def2tzvp
Making QCJSON for neb-bare.0-orca.sp-bp86.def2svp

0 file(s) encountered errors and not written
```

ORCA calculation examples descriptions.

- `4h2o.abc0.iter2-orca.opt-mp2.def2tzvp.out`: An optimization using MP2/def2-TZVP with the default frozen core approximation.
- `neb-bare.0-orca.sp-bp86.def2svp.out`: A single-point energy calculation using BP86-D3BJ/def2-TZVP.
- `neb-bare.0-orca.sp-bp86.def2tzvp.smd.out`: A single-point energy calculation using BP86-D3BJ/def2-TZVP with a water implicit solvent model (SMD).
- `neb-bare.0-orca.sp-wb97x.def2tzvp.out`: A single-point energy calculation using $\omega$B97X-D3BJ/def2-TZVP.
- `5h2o.abc0.iter1.mp2.md.300k.iter1.mol0,1,2,3,4-orca.engrad-mp2.def2tzvp.out`: An ORCA job running multiple energy+gradient calculations on different configurations of the same system.
- `alo4h5-orca.freq-bp86.def2svp.out`: An analytical frequency calculation with zero imaginary frequencies.
- `alo3h3.h2o-orca.freq-bp86.def2svp.cpcm.out`: An analytical frequency calculation with an implicit solvent and one imaginary frequency.

## QCJSON Format

Custom modifications are made to [QCSchema](https://github.com/MolSSI/QCSchema) to meet the immediate needs of our research.

- [QCSchema](https://github.com/MolSSI/QCSchema) does not support multiple configurations in a single JSON file: one file per structure and calculation.
  This is inherently incompatible with geometry optimizations, trajectories, or any other iterative procedure.
  At of the time of writing (2021-02-09), there has been no consensus of how to manage these files.
  Our immediate solution is to list each [QCSchema](https://github.com/MolSSI/QCSchema) inside a list (i.e., [ { }, { }, { }, ... ] ).
  Furthermore, only iterations that have the same topology are supported (no checks are provided).
- Other information with keywords not defined in the [QCSchema](https://github.com/MolSSI/QCSchema).

### Keys

Here are definitions of all [QCSchema](https://github.com/MolSSI/QCSchema) and custom keys used in this script.
Not all QCJSONs will have every key all the time; some are method dependent like the frozen core approximation.
All custom keys are marked with a *, and the Python type is specified.

- ``"schema_name"``: str, specifies the type of QCJSON. Always equal to "qc_schema_output".
- ``"schema_version"``: int, specifies the [QCSchema](https://github.com/MolSSI/QCSchema) version.
  Current version is 1.
- ``"qcjson_creator_version"``*: str, specifies the qcjson-creator script version.
- ``"provenance"``: dict, a brief description of the program, version, and routine used to generate the output.
  Can include more detailed information such as computation time, processor information, and host location.
    - ``"creator"``: str, name of the QC package.
    - ``"version"``: str, version of the QC package.
- ``"molecule"``: dict, a full description of the overall molecule(s) (e.g., its geometry, fragments, and charges).
    - ``"geometry"``: list [[float, float, float], [float, float, float], ...], vector of XYZ coordinates of the atoms of length equal to the number of atoms.
    - ``"symbols"``: list [str, str, ...], atom symbols in title case (e.g., "H" and "Na").
- ``"molecular_charge"``: int, the overall charge of the molecule.
- ``"molecular_multiplicity"``: int, the overall multiplicity of the molecule.
- ``"name"``: str, desired name of the molecule or system (currently is just the file name).
- ``"atomic_numbers"``: list [int, int, ...], atomic numbers, nuclear charge for atoms.
  Ghostedness should be indicated through ‘real’ field, not zeros here.
- ``"driver"``: str, what are you looking to calculate: energy, gradient, Hessian, or property.
  Note, we implement other driver options: optimization.
- ``"return_result"``: str, the “primary” return of a given computation.
  For energy, gradient, and Hessian quantities these are either single numbers or a array representing the derivative quantity.
- ``"model"``: dict, the overall mathematical model we are using for our calculation.
  Another way to think about this is the largest superset that still obtains roughly the same result.
  In QM, this is the Hamiltonian (HF, DFT, ...) combined with the overall basis of the calculation.
  An example in QM would be HF/STO-3G or B3LYP/6-311G**.
  Custom basis sets can be handled with custom keywords.
    - ``"method"``: str, primary method or level of theory.
    - ``"basis"``: str, overall basis set.
    - ``"aux_basis"``*: str, auxillary basis set if RI approximation is used.
- ``"keywords"``: dict, various tunable parameters for the calculation.
  These vary widely, depending on the basis and model chemistry.
  These represent the keywords of individual programs currently.
    - ``"scf_convergence_tolerance"``*: str, program-specific SCF convergence criteria.
      For example, "tight".
    - ``"dispersion"``*: str, if empirical dispersion corrections are used (with DFT), this specifies the method.
    - ``"scf_grid_level"``*: int, specifies program-specific integration grid level for the scf cycle.
    - ``"final_grid_level"``*: int, specifies program-specific integration grid level for the final energy calculation.
    - ``"implicit_solvent"``*: str, the implicit solvent model if used.
      For example, "SMD" or "CPCM".
      Note, this property is only included if included.
    - ``"solvent_name"``*: str, name of the solvent.
      For example, "water".
      Note, this property is only included if an implicit solvent is used.
    - ``"rij_approximation"``*: bool, if the resolution of identity approximation for the Coulomb integrals (RI-J) is used.
      Note, this property is only included if true.
    - ``"cosx_approximation"``*: bool, if the chain-of-spheres integration approximation to the exchange term (COSX) is used.
      Note, this property is only included if true.
- ``"properties"``: dict, a list of valid quantum chemistry properties tracked by the schema.
    - ``"calcinfo_nbasis"``: int, the number of basis functions for the computation.
    - ``"calcinfo_nmo"``: int, the number of molecular orbitals for the computation.
    - ``"return_energy"``: float, the energy of the requested method, identical to return_value for energy computations.
      For frequency calculations, this will be the Gibbs free energy.
    - ``"scf_total_energy"``: float, the total electronic energy of the SCF stage of the calculation.
      This is represented as the sum of the ... quantities.
    - ``"scf_iterations"``: int, the number of SCF iterations taken before convergence.
    - ``"scf_dispersion_correction_energy"``: float, the dispersion correction appended to an underlying functional when a DFT-D method is requested.
    - ``"scf_one_electron_energy"``: float, the one-electron (core Hamiltonian) energy contribution to the total SCF energy.
    - ``"scf_two_electron_energy"``: float, the two-electron energy contribution to the total SCF energy.
    - ``"scf_xc_energy"``: float, the functional energy contribution to the total SCF energy.
    - ``"scf_dipole_moment"``: list [float, float, ...], the x, y, and z dipole components.
    - ``"nuclear_repulsion_energy"``: float, the nuclear repulsion energy contribution to the total SCF energy.
    - ``"mp2_total_energy"``: float, the total MP2 energy (MP2 correlation energy + HF energy).
    - ``"mp2_correlation_energy"``: float, the MP2 correlation energy.
    - ``"mp2_dipole_moment"``: list [float, float, ...], the x, y, and z dipole components.
    - ``"alpha_homo_energy"``*: float, highest occupied molecular orbital energy of the alpha electron.
    - ``"alpha_homo_lumo_gap_energy"``*: float, energy difference between lowest unoccupied and highest occupied molecule orbital of the alpha electron.
    - ``"beta_homo_energy"``*: float, highest occupied molecular orbital energy of the beta electron.
    - ``"beta_homo_lumo_gap_energy"``*: float, energy difference between lowest unoccupied and highest occupied molecule orbital of the alpha electron.
    - ``"mulliken_charges"``*: list [float, float, ...], the Mulliken atomic charges.
    - ``"loewdin_charges"``*: list [float, float, ...], the Loewdin atomic charges.
    - ``"omega"``*: list [float, float, ...], vibrational frequencies in cm<sup>-1</sup>.
      Includes zero frequencies.
    - ``"q"``*: list [[float, float, ...], [float, float, ...], ...], the normalized, mass-weighted normal modes.
    - ``"temperature"``*: float, the temperature, in Kelvin, used for thermochemistry.
    - ``"zero_point_vibrational_correction"``*: float, the zero-point vibrational energy correction at 0 K to the electronic energy.
    - ``"thermal_energy_corrections"``*: float, the thermal vibrational, rotational, and translational correction to the electronic energy.
    - ``"enthalpic_corrections"``*: float, the thermal enthalpy corrections for the electronic energy at the specified temperature.
    - ``"entropic_corrections"``*: float, the translational, rotational, and vibrational entropic corrections to enthalpy to get Gibbs free energy (i.e., TS).
- ``"success"``: bool, a description if the computation was successful or not.
- ``"error"``: dict, for unsuccessful computations standard errors will be placed in the output such as convergence, IO errors, etc.
    - ``"error_type"``: str, what caused the error.
      Not yet implemented in this script.
    - ``"error_message"``: str, specific program message.
      Not yet implemented in this script.

## License

Distributed under the MIT License. See `LICENSE` for more information.
