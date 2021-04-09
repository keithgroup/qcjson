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
from qcjson.utils import atoms_by_number
from qcjson.utils import convert_forces
from qcjson.utils import parse_stringfile
from qcjson.jsons import QCJSON
from qcjson.parsers import xtbParser

class xtbJSON(QCJSON):
    """xtb specific QCJSON information.

    Supported calculation types:

    - Optimization
    - Energies
    - Molecular dynamics

    Parameters
    ----------
    outfile_path : :obj:`str`
        Path to output file.
    geomfile_path : :obj:`str`
        Path to file containing xyz structures of the calculation.

    Attributes
    ----------
    outfile_path : :obj:`str`
        Path to output file.
    method_type : :obj:`str`
        QCSchema method type of ``'scf'`` (e.g., DFT), ``'moller-plesset'``,
        or ``'coupled cluster'``. For xtb this is always ``'scf'``.
    parser : :obj:`xtbParser`
        Manually parse information from output file.
    """

    def __init__(self, outfile_path, geomfile_path):
        super().__init__()
        self.outfile_path = outfile_path
        self.geomfile_path = geomfile_path
        self.parse_output()
        self.method_type = None
    
    def parse_output(self):
        """Parse output file using cclib and custom parser.
        """
        filename_with_extension = os.path.basename(self.outfile_path)
        self.path = os.path.abspath(self.outfile_path)
        self.name = '.'.join(filename_with_extension.split('.')[:-1])
        
        # Custom parsed information.
        self.parser = xtbParser(self.path)
        self.parser.parse()
        self.parsed_data = self.parser.data

        # Gets geometries.
        z, comments, data = parse_stringfile(self.geomfile_path)
        self._symbols = z
        self._atomic_numbers = [atoms_by_number(i) for i in z]
        self._geometry = data

        # In molecular dynamics simulation with xtb, energy dumps are not
        # for every structure given in the trajectory.
        # The default energy parsing will parse energies for the first structure.
        # For optimizations, the energies have more significant figures in the
        # geometry file.
        # So we will overwrite our scf_total_energies with those parsed from
        # the geometry file.
        # We do this here and not in xtbParser to avoid parsing twice.
        if self.parsed_data['driver'] in ['optimization', 'molecular dynamics']:
            energies = []
            for comment in comments:
                split = comment.split()
                energies.append(float(split[1]))
            self.parsed_data['properties']['scf_total_energy'] = energies
    
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
                    'geometry' : self._geometry[iteration],
                    'symbols': self._symbols[iteration]
                },
                'molecular_charge': self.parsed_data['molecular_charge'],
                'molecular_multiplicity': self.parsed_data['molecular_multiplicity'],
                'name': self.name,
                'atomic_numbers': self._atomic_numbers[iteration]
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
        """
        model = {
            'method': self.parsed_data['model']['method']
        }
        
        return model
    
    def get_keywords(self, iteration=-1):
        """Package-specific job properties.

        This will be placed under the ``'keyword'`` JSON property.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        """
        keywords = {}

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
            Calculated properties of the system.
        """
        properties = {}

        # Some properties are not parsed for every structure and cause
        # errors with the iteration scheme we use.
        # For now we will not include these properties.

        # Nuclear repulsion energy is only printed for the first and
        # last structures during optimizations.
        opt_skip = ['nuclear_repulsion_energy']
        # Same issue for molecular dynamics.
        md_skip = ['nuclear_repulsion_energy']

        # Add parsed property information.
        for info in self.parsed_data['properties'].keys():
            if self.parsed_data['driver'] == 'optimization':
                if info in opt_skip:
                    continue
            
            if self.parsed_data['driver'] == 'molecular dynamics':
                if info in md_skip:
                    continue
            
            data = self.parsed_data['properties'][info][iteration]
            properties[info] = data
        
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
        driver = {
            'driver': self.parsed_data['driver']
        }
        
        # Adds return_result for supported calculations.
        if self.parsed_data['driver'] == 'energy':
            return_result = self.parsed_data['properties']['scf_total_energy'][iteration]
            driver['return_result'] = return_result
        
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
            'creator': 'xtb',
            'version': self.parsed_data['provenance']['version']
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
        # pylint: disable=undefined-variable
        if not hasattr(self, '_json'):

            all_jsons = []
            
            # Creates QCJSON for each geometry in the calculation.
            for i in range(0, len(self._geometry)):
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
                    all_jsons[-1]['success'] = self.parsed_data['success']
                    if len(all_jsons) == 1:
                        self._json = all_jsons[0]
                    else:
                        self._json = all_jsons
                except Exception:
                    if debug:
                        raise
                    else:
                        # These functions and variables are in the qcjson-creator
                        # script.
                        if self.path not in error_files:
                            error_out(self.path, 'Uncaught exceptions.')
                        self._json = all_jsons
                        break
        
        return self._json
