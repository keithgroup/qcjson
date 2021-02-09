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

import cclib
import json
from qcjson import __version__ as qcjson_version

class QCJSON:
    """Base quantum chemistry JSON class.

    Attributes
    ----------
    outfile : :obj:'cclib.parser.data.ccData_optdone_bool'
        Parsed data from output file.
    name : :obj:`str`
        Name identifying the QCJSON. Defaults to ``'qcschema'`` until an output
        file is parsed by cclib.
    path : :obj:`str`
        Absolute path to the output file.
    multiple_iterations : :obj:`bool`
        If the QCJSON is from a calculation with multiple iterations.
    """

    def __init__(self):
        self.name = 'qcjson'
        self.cclib_data = None
        self.path = None

    def write(self, name, json_dict, save_dir, prettify=True):
        """Writes JSON.

        Parameters
        ----------
        name : :obj:`str`
            Name of file.
        save_dir : :obj:`str`
            Path to save directory.
        prettify : :obj:`bool`
            Indents JSON objects if True. If false the JSON file is only one
            line.
        """
        if save_dir[-1] != '/':
            save_dir += '/'
        if prettify:
            json_string = json.dumps(
                json_dict, cls=cclib.io.cjsonwriter.JSONIndentEncoder,
                sort_keys=True, indent=4
            )
        else:
            json_string = json.dumps(
                json_dict, cls=cclib.io.cjsonwriter.NumpyAwareJSONEncoder,
                sort_keys=True
            )
        with open(f'{save_dir}{name}.json', 'w') as f:
            f.write(json_string)
    
    @property
    def schema(self):
        """Base QCSchema information.
        """
        if not hasattr(self, 'cclib_data'):
            raise AttributeError('No output file was parsed.')
        schema_dict = {
            'schema_name': 'qc_schema_output',
            'schema_version': 1,
            'qcjson_creator_version': qcjson_version
        }
        return schema_dict
    
    def _get_scf_energy(self, iteration=-1):
        """The energy after SCF cycle.

        This is the energy for all DFT methods and HF energy for post-HF
        calculations.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`float`
            SCF energy in Hartree.
        """
        scfenergy = cclib.parser.utils.convertor(
            self.cclib_data.scfenergies[iteration], 'eV', 'hartree'
        )
        return scfenergy
    
    def _get_dispersion_energy(self, iteration=-1):
        """Dispersion energy corrections to DFT calculations.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`float`
            Dispersion energy correction in Hartree.
        """
        dispersion_energy = cclib.parser.utils.convertor(
            self.cclib_data.dispersionenergies[iteration], 'eV', 'hartree'
        )
        return dispersion_energy
    
    def _get_mp_energy(self, iteration=-1):
        """Total energies after Moller-Plesset corrections.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`float`
            Total energy from highest order MP corrections.
        """
        if self.cclib_data.mpenergies.ndim == 1:
            mpenergy_ev = self.cclib_data.mpenergies[iteration]
        elif self.cclib_data.mpenergies.ndim == 2:
            order = self.cclib_data.mpenergies.shape[1] + 1
            mpenergy_ev = self.cclib_data.mpenergies[iteration][order - 2]
        mpenergy_hartree = cclib.parser.utils.convertor(
            mpenergy_ev, 'eV', 'hartree'
        )
        return mpenergy_hartree
    
    def _get_cc_energy(self, iteration=-1):
        """Energy after all coupled cluster corrections.

        Only the highest order correction is included.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`float`
            CC energy in Hartree.
        """
        ccenergy = cclib.parser.utils.convertor(
            self.cclib_data.ccenergies[iteration], 'eV', 'hartree'
        )
        return ccenergy

