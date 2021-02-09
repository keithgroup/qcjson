# MIT License
# 
# Copyright (c) 2021, Alex M. Maldonado
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import cclib
from packaging import version
from qcjson.utils import methods
from qcjson.utils import basis_sets
from qcjson.utils import atoms_by_element
from qcjson.utils import convert_forces
from qcjson.jsons.qcjson import QCJSON
from qcjson.parsers.orcaparser import orcaParser

class orcaJSON(QCJSON):
    """ORCA specific QCJSON information.

    Supported calculation types:

    - Single-point energies
    - Energy+Gradient
    - Optimizations

    Parameters
    ----------
    outfile_path : :obj:`str`
        Path to output file.

    Attributes
    ----------
    outfile_path : :obj:`str`
        Path to output file.
    method_type : :obj:`str`
        QCSchema method type of ``'scf'`` (e.g., DFT), ``'moller-plesset'``,
        or ``'coupled cluster'``.
    parser : :obj:`orcaParser`
        Manually parse information from output file.
    """

    def __init__(self, outfile_path):
        super().__init__()
        self.cclib_data = None
        self.outfile_path = outfile_path
        self.parse_output()
    
    def parse_output(self):
        """Parse output file using cclib and custom parser.
        Parameters
        ----------
        outfile_path : :obj:`str`
            Path to computational chemistry output file.
        """
        filename_with_extension = os.path.basename(self.outfile_path)
        self.path = os.path.abspath(self.outfile_path)
        self.name = '.'.join(filename_with_extension.split('.')[:-1])

        # cclib parsed information.
        self.cclib_data = cclib.io.ccread(self.outfile_path)
        if self.cclib_data.atomcoords.shape[0] > 1:
            self.multiple = True
        elif self.cclib_data.atomcoords.shape[0] == 1:
            self.multiple = False
        
        # Custom parsed information.
        self.parser = orcaParser(self.path)
        self.parser.parse()
        self.parsed_data = self.parser.data
    
    def get_topology(self, iteration=-1):
        """A full description of the overall molecule its geometry, fragments,
        and charges.

        Returned keys will be a top-level JSON property.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`dict`

            ``'molecule'``
                Atomic information about the system. Contains the following keys

                ``'geometry'``
                    (3 - nat, ) vector of XYZ coordinates [a0] of the atoms.
                ``'symbols'``
                    (nat, ) atom symbols in title case.
            ``'molecular_charge'``
                The overall charge of the molecule.
            ``'molecular_multiplicity'``
                The overall multiplicity of the molecule.
            ``'name'``
                The name of the molecule or system.
            ``'atomic_numbers'``
                (nat, ) atomic numbers, nuclear charge for atoms. Ghostedness
                should be indicated through ‘real’ field, not zeros here.
        """
        try:
            topology = {
                'molecule': {
                    'geometry' : self.cclib_data.atomcoords[iteration],
                    'symbols': atoms_by_element(self.cclib_data.atomnos.tolist())
                },
                'molecular_charge': self.cclib_data.charge,
                'molecular_multiplicity': self.cclib_data.mult,
                'name': self.name,
                'atomic_numbers': self.cclib_data.atomnos
            }
        except:
            raise ValueError('Topology data not succussfully parsed.')
        return topology

    def get_model(self, iteration=-1):
        """Model chemistry used for the calculation.

        This will be placed under the ``'model'`` JSON property.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`dict`

            ``'method'``
                They main quantum-chemical method used for the calculation. For
                DFT calculations this would be the functional.
            ``'basis'``
                The basis set used for the calculation.
        """
        model = {}
        _remove_keywords = []
        for kw in self.orca_keywords:
            kw_lower = kw.lower()
            if kw_lower in methods['orca']['scf']:
                self.method_type = 'scf'
                # Ensures dispersion method is included in 'properties' for
                # functionals that automatically include it.
                model['method'] = kw
                if kw_lower == 'wb97x-d3' or kw_lower == 'b2plyp-d3':
                    self.orca_keywords.append('d3')
                if  kw_lower == 'b97m-d3bj' or kw_lower == 'wb97x-d3bj':
                    self.orca_keywords.append('d3bj')
                _remove_keywords.append(kw)
                break
            elif kw_lower in methods['orca']['moller-plesset']:
                self.method_type = 'moller-plesset'
                if 'mp3' in kw_lower:
                    raise ValueError('MP3 is not supported.')
                model['method'] = kw
                _remove_keywords.append(kw)
                break
            elif kw_lower in methods['orca']['coupled cluster']:
                self.method_type = 'coupled cluster'
                model['method'] = kw
                _remove_keywords.append(kw)
                break
        
        # Need a separate keyword iteration for basis sets because
        # dispersion-corrected functionals messes up the indicies for pop.
        for kw in self.orca_keywords:
            kw_lower = kw.lower()
            if kw_lower in basis_sets['orca']:
                model['basis'] = kw
                _remove_keywords.append(kw)
                break
        
        # Handles auxiliary basis sets parsed from output file.
        # Removes explicit keyword for specifying aux basis sets (will be
        # included from manually parsed custom information).
        for kw in self.orca_keywords:
            kw_lower = kw.lower()
            if kw_lower in ['def2/j']:
                _remove_keywords.append(kw)
        
        # Adds manually parsed model information.
        for info in self.parsed_data['model'].keys():
            data = self.parsed_data['model'][info]
            if type(data) == list:
                if len(data) != 1:
                    # ORCA does an additional evaluation that does not have
                    # geometric convergence information.
                    if len(data) == iteration:
                        continue
                    else:
                        data = data[iteration]
            model = {
                **model, **{info: data}
            }
        
        for kw in _remove_keywords:
            self.orca_keywords.remove(kw)
        
        return model
    
    def get_keywords(self, iteration=-1):
        """Package-specific job properties.

        This will be placed under the ``'keyword'`` JSON property.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`dict`
            Could contain the following keys:

            ``'dispersion'``
                Method of empirical dispersion included in job.
            ``'implicit_solvent'``
                Implicit (i.e., continuum) solvent model used in job.
            ``'solvent_name'``
                Name of the solvent (e.g., ``'water'``).
            ``'frozen_core'``
                If the FrozenCore approximation was used. This defaults to
                ``True`` for MP and CC calculations in ORCA.
            ``'scf_convergence_tolerance'``
                Self-consistent field convergence tolerance keyword for ORCA.
            ``'rij_approximation'``
                The resolution of identity (RI) approximation for the Coulomb
                (J) term.
            ``'cosx_approximation'``
                The chain-of-spheres integration approximation to the exchange
                term (COSX).
            ``'scf_grid_level'``
                :obj:`int` specifying ORCA default grid level (1 to 7) during
                the SCF cycle.
            ``'final_grid_level'``
                :obj:`int` specifying ORCA default grid level for the final
                evaluation.
        """
        keywords = {}
        _remove_keywords = []
        
        # Parses and categorizes calculation parameters.
        for kw in self.orca_keywords:
            kw_lower = kw.lower()

            # Empirical dispersion
            if kw_lower in ['d4', 'd3bj', 'd3', 'd3zero', 'd2']:
                if version.parse(self.cclib_data.metadata['package_version']) \
                   >= version.parse('4.0.0') and kw_lower == 'd3':
                    # In ORCA 4, D3 is D3BJ
                    keywords['dispersion'] = 'D3BJ'
                else:
                    keywords['dispersion'] = kw
                if kw not in _remove_keywords:
                    _remove_keywords.append(kw)
                else:
                    break
            
            # Implicit solvent models
            if 'cpcm' in kw_lower or 'c-pcm' in kw_lower:
                if 'smd' in self.cclib_data.metadata['input_file_contents']:
                    keywords['implicit_solvent'] = 'SMD'
                else:
                    keywords['implicit_solvent'] = 'CPCM'
                if '(' in kw_lower and ')' == kw_lower[-1]:
                    solvent_name = kw_lower[:-1].split('(')[-1]
                    keywords['solvent_name'] = solvent_name
                if kw not in _remove_keywords:
                    _remove_keywords.append(kw)
                else:
                    break
            
            if kw_lower == 'frozencore':
                keywords['frozencore'] = True
                if kw not in _remove_keywords:
                    _remove_keywords.append(kw)
                else:
                    break
            elif kw_lower == 'nofrozencore':
                keywords['frozencore'] = False
                if kw not in _remove_keywords:
                    _remove_keywords.append(kw)
                else:
                    break
            
            if kw_lower in ['sloppyscf', 'loosescf', 'normalscf', 'strongscf',
                            'tightscf', 'verytightscf', 'extremescf'] \
                or 'scfconv' in kw_lower:
                keywords['scf_convergence_tolerance'] = kw[:-3]
                if kw not in _remove_keywords:
                    _remove_keywords.append(kw)
                else:
                    break
            
            # Removes SCF information (will be manually parsed from outfile).
            if 'ri' == kw_lower or 'rijcosx' in kw_lower or 'rijk' == kw_lower \
               or 'ri-jk' == kw_lower:
                if kw not in _remove_keywords:
                    _remove_keywords.append(kw)
                else:
                    break
            
            # Removes grid information (will be manually parsed from outfile).
            if 'grid' in kw_lower:
                if kw not in _remove_keywords:
                    _remove_keywords.append(kw)
                else:
                    break
        
        for kw in _remove_keywords:
            self.orca_keywords.remove(kw)
        
        # Add default scf convergence if not already specified.
        if 'scf_convergence_tolerance' not in keywords.keys():
            keywords['scf_convergence_tolerance'] = 'Normal'
        
        # Adds manually parsed information.
        for info in self.parsed_data['keywords'].keys():
            data = self.parsed_data['keywords'][info]
            if type(data) == list:
                if len(data) != 1:
                    # ORCA does an additional evaluation that does not have
                    # geometric convergence information.
                    if len(data) == iteration:
                        continue
                    else:
                        data = data[iteration]
            keywords = {
                **keywords, **{info: data}
            }
        
        # Specifies if FrozenCore is used in MP or CC jobs.
        if self.method_type == 'moller-plesset':
            if version.parse(self.cclib_data.metadata['package_version']) \
                >= version.parse('4.0.0'):
                lower_list = [i.lower() for i in self.orca_keywords]
                if 'nofrozencore' in lower_list:
                    keywords['frozen_core'] = False
                    self.orca_keywords.pop(lower_list.index('nofrozencore'))
                else:
                    if 'frozencore' not in keywords.keys():
                        keywords['frozen_core'] = True

        # Add uncategorized calculation properties
        if len(self.orca_keywords) != 0:
            keywords['other'] = self.orca_keywords

        return keywords
    
    def get_properties(self, iteration=-1):
        """A list of valid quantum chemistry properties tracked by the schema.

        This will be placed under the ``'properties'`` JSON property.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`dict`

            ``'calcinfo_nbasis'``
                The number of basis functions for the computation.
            ``'calcinfo_nmo'``
                The number of molecular orbitals for the computation.
            ``'scf_total_energy'``
                The total electronic energy of the SCF stage of the calculation.
                This is represented as the sum of the … quantities.
            ``'scf_dispersion_correction_energy'``
                The dispersion correction appended to an underlying functional
                when a DFT-D method is requested.
            ``'scf_iterations'``
                The number of SCF iterations taken before convergence.
            ``'mp2_total_energy'``
                The total MP2 energy (MP2 correlation energy + HF energy).
        """
        properties = {}
        if self.method_type == 'scf':
            properties['scf_total_energy'] = self._get_scf_energy(iteration=iteration)
            electronic_energy = properties['scf_total_energy']
            properties['scf_dispersion_correction_energy'] = cclib.parser.utils.convertor(
                self.cclib_data.dispersionenergies[iteration], 'eV', 'hartree'
            )
            properties['scf_iterations'] = self.cclib_data.scfvalues[0].shape[0]
            if 'dipole_moment' in self.parsed_data['properties'].keys():
                properties['scf_dipole_moment'] = self.parsed_data['properties']['dipole_moment']
        elif self.method_type == 'moller-plesset':
            properties['scf_total_energy'] = self._get_scf_energy(iteration=iteration)
            properties['mp2_total_energy'] = self._get_mp_energy(iteration=iteration)
            electronic_energy = properties['mp2_total_energy']

            # MP2 correlation energies are included in every iteration of jobs.
            # Thus, even with optimizations, we will include other mp2 energies.
            for info in self.parsed_data['properties'].keys():
                if 'mp2_' in info:
                    data = self.parsed_data['properties'][info][iteration]
                    data_merge = {info: data}
                    properties = {**properties, **data_merge}
        elif self.method_type == 'coupled cluster':
            # TODO
            pass
        else:
            raise ValueError('Unknown method type.')
        
        # Manages if return_energy is electronic energy or Gibbs free energy.
        if self.calc_driver == 'frequency':
            gibbs_corr = 0.0
            for correction in [
                'zero_point_vibrational_correction', 'thermal_energy_corrections',
                'enthalpic_corrections', 'entropic_corrections'
            ]:
                if correction == 'entropic_corrections':
                    gibbs_corr -= self.parsed_data['properties'][correction][iteration]
                else:
                    gibbs_corr += self.parsed_data['properties'][correction][iteration]
            properties['return_energy'] = electronic_energy + gibbs_corr
        else:
            properties['return_energy'] = electronic_energy

        # Other calculation information.
        properties['calcinfo_nbasis'] = self.cclib_data.nbasis
        properties['calcinfo_nmo'] = self.cclib_data.nmo

        # Adds manually parsed data.
        # ORCA optimizations does not always include all the components for
        # each iteration. So, we don't include the extra information here.
        if self.calc_driver != 'optimization':
            for info in self.parsed_data['properties'].keys():
                data = self.parsed_data['properties'][info][iteration]

                # Rename any properties
                if info == 'dipole_moment':
                    if self.method_type == 'scf':
                        prefix = 'scf'
                    elif self.method_type == 'moller-plesset':
                        prefix = 'mp2'
                    elif self.method_type == 'coupled cluster':
                        pass
                    info = prefix + '_' + info

                data_merge = {info: data}
                properties = {**properties, **data_merge}
        
        # Alpha and beta electron properties
        homo_idx_alpha = int(self.cclib_data.homos[0])
        homo_idx_beta = int(self.cclib_data.homos[-1])
        alpha_homo_energy = cclib.parser.utils.convertor(
            self.cclib_data.moenergies[0][homo_idx_alpha], 'eV', 'hartree'
        )
        alpha_lumo_energy = cclib.parser.utils.convertor(
            self.cclib_data.moenergies[0][homo_idx_alpha + 1], 'eV', 'hartree'
        )
        alpha_gap_energy = alpha_lumo_energy - alpha_homo_energy
        beta_homo_energy = cclib.parser.utils.convertor(
            self.cclib_data.moenergies[-1][homo_idx_beta], 'eV', 'hartree'
        )
        beta_lumo_energy = cclib.parser.utils.convertor(
            self.cclib_data.moenergies[-1][homo_idx_beta + 1], 'eV', 'hartree'
        )
        beta_gap_energy = beta_lumo_energy - beta_homo_energy
        properties['alpha_homo_energy'] = alpha_homo_energy
        properties['alpha_homo_lumo_gap_energy'] = alpha_gap_energy
        properties['beta_homo_energy'] = beta_homo_energy
        properties['beta_homo_lumo_gap_energy'] = beta_gap_energy
        
        return properties
    
    def get_driver(self, iteration=-1):
        """The purpose of the calculation and its direct result.

        Returned keys will be a top-level JSON property.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`dict`

            ``'driver'``
                The purpose of the calculation. Could be ``'energy'``,
                ``'gradient'``, ``'optimization'``, or ``'frequency'``.
            ``'return_results'``
                The direct result of the driver calculation. For example, a
                energy+gradient calculation will return the gradient (the
                energies will be included in the ``'properties'`` property).
        """
        driver = {}
        _remove_keywords = []

        for kw in self.orca_keywords:
            kw_lower = kw.lower()

            # Gradients
            if kw_lower == 'engrad' or kw_lower == 'numgrad':
                self.calc_driver = 'gradient'
                driver['driver'] = self.calc_driver

                if self.cclib_data.grads.ndim == 3:
                    grads = self.cclib_data.grads[iteration]
                else:
                    raise ValueError('Please check gradient dimensions.')
                driver['return_result'] = convert_forces(
                    grads, 'hartree', 'bohr', 'hartree', 'Angstrom'
                )

                _remove_keywords.append(kw)
                break

            # Frequencies (not supported)
            elif kw_lower == 'freq' or kw_lower == 'numfreq':
                self.calc_driver = 'frequency'
                driver['driver'] = self.calc_driver

                _remove_keywords.append(kw)
                break

            # Optimizations
            elif kw_lower == 'opt' or kw_lower == 'copt' or kw_lower == 'zopt':
                self.calc_driver = 'optimization'
                driver['driver'] = self.calc_driver

                _remove_keywords.append(kw)
                break
            
            # Energies
            elif kw_lower == 'energy' or kw_lower == 'sp':
                # This driver will be captured by the len(driver) == 0 condition.
                self.calc_driver = 'energy'
                _remove_keywords.append(kw)
                break
        
        # ORCA defaults to single-point energies if no keyword is present.
        if not hasattr(self, 'calc_driver') or self.calc_driver == 'energy':
            self.calc_driver = 'energy'
            driver['driver'] = self.calc_driver
            if hasattr(self.cclib_data, 'ccenergies'):
                return_result = self.cclib_data.ccenergies[iteration]
            elif hasattr(self.cclib_data, 'mpenergies'):
                return_result = self._get_mp_energy(
                iteration=iteration
            )
            elif hasattr(self.cclib_data, 'scfenergies'):
                return_result = self._get_scf_energy(
                    iteration=iteration
                )
            driver['return_result'] = return_result
            pass

        # Removes triggered keywords so as not to be included in uncategorized
        # keywords.
        for kw in _remove_keywords:
            self.orca_keywords.remove(kw)
        
        return driver

    def get_provenance(self):
        """Information about the originating package of the calculation.

        Returns
        -------
        :obj:`dict`
            
            ``'creator'``
                Package name.
            ``'version'``
                Package version.
        """
        provenance = {
            'creator': 'ORCA',
            'version': self.cclib_data.metadata['package_version']
        }
        return provenance
    
    def get_json(self, debug=False):
        """QCJSON of an ORCA output file.

        Calculations supported:
        - single-point energies
        - gradients
        - optimizations

        :type: :obj:`dict`
        """
        if not hasattr(self, '_json'):

            all_jsons = []
            self.orca_keywords = self.cclib_data.metadata['keywords'].copy()
            for i in range(0, self.cclib_data.atomcoords.shape[0]):
                # Optimizations are iterative with only a single set of
                # keywords; this is different than consecutive jobs (e.g., 
                # energies) that have repeated keywords. So, we reinitialzie the
                # keywords for each optimization iteration.
                if hasattr(self, 'calc_driver') \
                   and self.calc_driver == 'optimization':
                    self.orca_keywords = self.cclib_data.metadata['keywords'].copy()
                try:
                    all_jsons.append(super().schema)
                    all_jsons[-1]['provenance'] = self.get_provenance()
                    all_jsons[-1] = {
                        **all_jsons[-1], **self.get_topology(iteration=i)
                    }
                    all_jsons[-1] = {
                        **all_jsons[-1], **self.get_driver(iteration=i)
                    }
                    all_jsons[-1]['model'] = self.get_model(iteration=i)
                    all_jsons[-1]['keywords'] = self.get_keywords(iteration=i)
                    all_jsons[-1]['properties'] = self.get_properties(iteration=i)
                    all_jsons[-1]['success'] = self.cclib_data.metadata['success']
                    if len(all_jsons) == 1:
                        self._json = all_jsons[0]
                    else:
                        self._json = all_jsons
                except Exception:
                    if debug:
                        raise
                    else:
                        if self.path not in error_files:
                            error_out(self.path, 'Uncaught exceptions.')
                        self._json = all_jsons
                        break
        
        return self._json
