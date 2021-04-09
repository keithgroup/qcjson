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
from qcjson.parsers.parser import outfileParser
from qcjson.utils import parse_stringfile

class xtbParser(outfileParser):
    """Custom parser for xtb output files.
    """

    def __init__(self, outfile_path):
        super().__init__(outfile_path)
        self.multiple = True  # xtb calculations are always geom_split.
    
    def extract(self, outfile, line):
        """Extracts all possible information from trajectory file.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        #    * xtb version 6.3.3 (5b13467) compiled by 'ehlert@majestix' on 2020-09-17
        if '* xtb version' in line:
            self._extract_version(outfile, line)
        
        # program call               : xtb 5h2o.example.xyz.xyz --scc --gfn 2 --charge 0
        if 'program call               :' in line.strip():
            self._extract_run_type(outfile, line)
        
        # charge                     :                     0
        if 'charge                     :' in line.strip():
            self._extract_charge(outfile, line)
        
        # spin                       :                   0.0
        if 'spin                       :' in line.strip():
            self._extract_multiplicity(outfile, line)

        # :  Hamiltonian                  GFN2-xTB          :
        if ':  Hamiltonian ' in line.strip():
            self._extract_hamiltonian(outfile, line)

        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        # ::                     SUMMARY                     ::
        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        # :: total energy             -25.404338997931 Eh    ::
        # :: gradient norm              0.048763527905 Eh/a0 ::
        # :: HOMO-LUMO gap             12.772900340050 eV    ::
        # ::.................................................::
        # :: SCC energy               -25.580127846040 Eh    ::
        # :: -> isotropic ES            0.153316435656 Eh    ::
        # :: -> anisotropic ES         -0.020934020006 Eh    ::
        # :: -> anisotropic XC         -0.003785806991 Eh    ::
        # :: -> dispersion             -0.004190013801 Eh    ::
        # :: repulsion energy           0.175746000529 Eh    ::
        # :: add. restraining           0.000000000000 Eh    ::
        # :: total charge              -0.000000000000 e     ::
        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        if '::                     SUMMARY                     ::' in line.strip():
            self._extract_summary_energies(outfile, line)
        
        #  ----------------------------------------------------------- 
        # |                   =====================                   |
        # |                        A N C O P T                        |
        # |                   =====================                   |
        # |               Approximate Normal Coordinate               |
        # |                Rational Function Optimizer                |
        #  ----------------------------------------------------------- 

        #  ----------------------------------------------------------- 
        # |                       L-ANC optimizer                     |
        #  ----------------------------------------------------------- 
        if 'A N C O P T' == line.strip('| \n') or 'L-ANC optimizer' == line.strip('| \n'):
            self._extract_opt_data(outfile, line)
        
        if 'normal termination of xtb' == line.strip():
            self.data['success'] = True

    def _after_parse(self):
        """Checks to perform after parsing output file.
        """
        # During optimizations, xtb prints the last energy twice.
        # The last printed energy has more significant figures, so we will
        # get rid of the second to last one.
        # It is also unclear why the structure considered "CYCLE   1" is
        # missing from the geometry log. It goes from the provided structure
        # to the "CYCLE   2" structure. So we will also remove this energy.
        # Note that we overwrite this data later, but it is better to ensure
        # consistency than punting to later in the code.
        if self.data['driver'] == 'optimization':
            del self.data['properties']['scf_total_energy'][1]
            del self.data['properties']['scf_total_energy'][-2]
        
        # Sets success to fail if not true.
        if 'success' not in self.data.keys():
            self.data['success'] = False

    def _extract_version(self, outfile, line):
        """Version of xtb.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        line_split = line.split()
        version = line_split[3]
        self.data['provenance'] = {
            'version': version
        }
        line = next(outfile)
    
    def _extract_run_type(self, outfile, line):
        """Calculation type (e.g., opt, sp, grad).

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        if '--scc' in line:
            driver = 'energy'
        elif '--opt' in line:
            driver = 'optimization'
        elif '--grad' in line:
            driver = 'gradient'
        elif '--ohess' in line or '--hess' in line:
            driver = 'frequency'
        elif '--omd' in line or '--md' in line:
            driver = 'molecular dynamics'
        self.data['driver'] = driver
        line = next(outfile)

    def _extract_charge(self, outfile, line):
        """Overall system charge.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        _, _, charge = line.strip(' :').split()
        self.data['molecular_charge'] = int(charge)
    
    def _extract_multiplicity(self, outfile, line):
        """Overall system charge.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        _, _, spin = line.strip().split()
        multiplicity = 2 * float(spin) + 1
        self.data['molecular_multiplicity'] = int(multiplicity)
    
    def _extract_hamiltonian(self, outfile, line):
        """The xTB hamiltonian used for the calculation.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        _, hamiltonian, _ = line.strip(' :').split()
        self.data['model']['method'] = hamiltonian
        line = next(outfile)
    
    def _extract_summary_energies(self, outfile, line):
        """Extracts energies listed in SUMMARY box.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        if 'scf_total_energy' not in self.data['properties'].keys():
            self.data['properties']['scf_total_energy'] = []
        
        if 'nuclear_repulsion_energy' not in self.data['properties'].keys():
            self.data['properties']['nuclear_repulsion_energy'] = []
        
        while '' != line.strip():
            if 'total energy' in line:
                line_split = line.split()
                scf_total_energy = float(line_split[3])
            
            if 'repulsion energy' in line:
                line_split = line.split()
                nuclear_repulsion_energy = float(line_split[3])

            line = next(outfile)

        self.data['properties']['scf_total_energy'].append(
            scf_total_energy
        )
        self.data['properties']['nuclear_repulsion_energy'].append(
            nuclear_repulsion_energy
        )

    def _extract_opt_data(self, outfile, line):
        """All information during optimization routine. Incldues setup, energy.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        
        # The following lines signal the termination of the optimization routine.
        # *** GEOMETRY OPTIMIZATION CONVERGED AFTER 73 ITERATIONS ***
        # *** FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN 200 CYCLES ***
        # Note that some lines have '(*******%)' in them, so we add spaces around
        # The asterisks to avoid triggering early.
        while ' *** ' not in line and 'GEOMETRY OPTIMIZATION' not in line:
            if '* total energy  :' in line:
                line_split = line.split()
                scf_total_energy = float(line_split[4])
                self.data['properties']['scf_total_energy'].append(
                    scf_total_energy
                )
            line = next(outfile)
