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
from qcjson.parsers.parser import outfileParser

class orcaParser(outfileParser):
    """Custom parser for ORCA output files.
    """

    def __init__(self, outfile_path):
        super().__init__(outfile_path)
    
    def extract(self, outfile, line):
        """Extracts all possible information from trajectory file.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        # ----- AuxJ basis set information -----
        # Your calculation utilizes the auxiliary basis: def2/J
        if 'Your calculation utilizes the auxiliary basis: ' in line.strip():
            self._extract_aux_basis(outfile, line)

        # -------------------
        # DFT GRID GENERATION
        # -------------------
        if 'DFT GRID GENERATION' == line.strip():
            self._extract_grid_info(outfile, line)
        
        # Setting up the final grid:
        if 'Setting up the final grid:' == line.strip():
            self._extract_grid_info(outfile, line)
        
        # ------------
        # SCF SETTINGS
        # ------------
        if 'SCF SETTINGS' == line.strip():
            self._extract_scf_info(outfile, line)

        # ----------------
        # TOTAL SCF ENERGY
        # ----------------
        if 'TOTAL SCF ENERGY' == line.strip():
            self._extract_scf_energies(outfile, line)
        
        #                       .--------------------.
        # ----------------------|Geometry convergence|-------------------------
        # Item                value                   Tolerance       Converged
        # ---------------------------------------------------------------------
        # Energy change      -0.0000629945            0.0000050000      NO
        # RMS gradient        0.0000464055            0.0001000000      YES
        # MAX gradient        0.0001724583            0.0003000000      YES
        # RMS step            0.0061958646            0.0020000000      NO
        # MAX step            0.0240077373            0.0040000000      NO
        # ........................................................
        # Max(Bonds)      0.0127      Max(Angles)    0.30
        # Max(Dihed)        0.32      Max(Improp)    0.00
        # ---------------------------------------------------------------------

        if 'Geometry convergence' in line and '----------------------' in line:
            self._extract_geo_conv(outfile, line)
        
        # ----------------------------------------------------------
        #                         ORCA  MP2 
        # ----------------------------------------------------------
        if 'ORCA  MP2' == line.strip():
            self._extract_mp_energies(outfile, line)
        
        # -----------------------
        # MULLIKEN ATOMIC CHARGES
        # -----------------------
        # 0 B :   -0.476812
        # 1 C :    0.422752
        # 2 H :   -0.129344
        # 3 H :   -0.117864
        if 'MULLIKEN ATOMIC CHARGES' == line.strip():
            self._extract_mulliken_charges(outfile, line)
        
        # ----------------------
        # LOEWDIN ATOMIC CHARGES
        # ----------------------
        # 0 B :   -0.996520
        # 1 C :   -0.397949
        # 2 H :    0.016421
        # 3 H :    0.006312      
        if 'LOEWDIN ATOMIC CHARGES' == line.strip():
            self._extract_loewdin_charges(outfile, line)
        
        # -------------
        # DIPOLE MOMENT
        # -------------
        #                                 X             Y             Z
        # Electronic contribution:      5.09504      -0.27709     -14.55675
        # Nuclear contribution   :     -3.13135       0.80269      11.24356
        #                         -----------------------------------------
        # Total Dipole Moment    :      1.96369       0.52560      -3.31319
        #                         -----------------------------------------
        # Magnitude (a.u.)       :      3.88710
        # Magnitude (Debye)      :      9.88023
        if 'DIPOLE MOMENT' == line.strip():
            self._extract_dipole(outfile, line)
        
        # -----------------------
        # VIBRATIONAL FREQUENCIES
        # -----------------------

        # Scaling factor for frequencies =  1.000000000  (already applied!)

        # 0:         0.00 cm**-1
        # 1:         0.00 cm**-1
        # 2:         0.00 cm**-1
        if line[:23] == 'VIBRATIONAL FREQUENCIES':
            self._extract_frequencies(outfile, line)

        # ------------
        # NORMAL MODES
        # ------------

        # These modes are the Cartesian displacements weighted by the diagonal matrix
        # M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom
        # Thus, these vectors are normalized but *not* orthogonal

        #                 0          1          2          3          4          5    
        #     0       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
        #     1       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
        if line[:12] == "NORMAL MODES":
            self._extract_normal_modes(outfile, line)
        
        # Temperature         ... 298.15 K
        if 'Temperature         ...' in line:
            self._extract_temp(outfile, line)
        
        # Zero point energy                ...      0.06455884 Eh      40.51 kcal/mol
        if 'Zero point energy                ...' in line:
            self._extract_zpve(outfile, line)
        
        # Total thermal correction                  0.00868494 Eh       5.45 kcal/mol
        if 'Total thermal correction' in line:
            self._extract_thermal_corrections(outfile, line)
        
        # Thermal Enthalpy correction       ...      0.00094421 Eh       0.59 kcal/mol
        if 'Thermal Enthalpy correction       ...' in line:
            self._extract_enthalpic_corr(outfile, line)
        
        # Final entropy term                ...      0.04136711
        if 'Final entropy term                ...' in line:
            self._extract_entropic_corr(outfile, line)

    def _after_parse(self):
        """Checks to perform after parsing output file.
        """

        if 'final_grid_level' not in self.data['keywords'].keys():
            self.data['keywords']['final_grid_level'] = self.data['keywords']['scf_grid_level']
        
        if 'geo_energy_change_value' in self.data['keywords'].keys():
            self.data['keywords']['geo_energy_change_value'].insert(0, 0.0)
    
    def _extract_aux_basis(self, outfile, line):
        """Information about auxiliary basis sets used in the calculation.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        line_split = line.strip().split(':')
        self.data['model']['aux_basis'] = line_split[-1].strip()


    def _extract_grid_info(self, outfile, line):
        """DFT integration grid information.

        All DFT calculations are performed with numerical integration, which
        means a integration grid must be specified. This is typically handeled
        with the Grid keyword. Furthermore, ORCA defaults to a multigrid
        approach where one grid is used for the SCF cycle, and a different
        (typically larger) grid is used for the final energy evaluation.

        After testing different keyword combinations, the Grid keywords are not
        always consistent with the final grid. Thus, we are going to directly
        parse the grid information from the output file instead of depending
        on the keywords.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        
        Notes
        -----
        Here are the main parameters specifying the ORCA default grid levels.
        Default SCF grid is 2.

        +--------+---------------+----------+
        |  Grid  |  AngularGrid  |  IntAcc  |
        +========+===============+==========+
        |  1     |  Lebedev-50   |  4.34    |
        |  2     |  Lebedev-110  |  4.34    |
        |  3     |  Lebedev-194  |  4.34    |
        |  4     |  Lebedev-302  |  4.67    |
        |  5     |  Lebedev-434  |  5.01    |
        |  6     |  Lebedev-590  |  5.34    |
        |  7     |  Lebedev-770  |  5.67    |
        +--------+---------------+----------+
        """
        lebedev_to_level = {
            '50': 1, '110': 2, '194': 3, '302': 4, '434': 5, '590': 6, '770': 7
        }
        
        # SCF Cycle Grid
        if 'DFT GRID GENERATION' == line.strip():
            while 'Angular Grid (max. acc.)' not in line:
                line = next(outfile)
            lebedev_num = line.strip().split('-')[-1]
            self.data['keywords']['scf_grid_level'] = lebedev_to_level[lebedev_num]
        elif 'Setting up the final grid:' == line.strip():
            while 'Angular Grid (max. acc.)' not in line:
                line = next(outfile)
            lebedev_num = line.strip().split('-')[-1]
            self.data['keywords']['final_grid_level'] = lebedev_to_level[lebedev_num]
    
    def _extract_scf_energies(self, outfile, line):
        """The nulear repulsion, one- and two-electron energy, and
        exchange-correlation energy after a SCF cycle.

        This is called directly after the ``'TOTAL SCF ENERGY'`` trigger, and 
        will terminate once the ``'SCF CONVERGENCE'`` trigger is reached.

        Instead of returning the energies themselves, we handle the creation and
        modification of ``scf_info`` here so any missing information (such as
        ``'scf_xc_energy'`` in MP2 calculations) is not an issue.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        while 'SCF CONVERGENCE' != line.strip():
            # Nuclear Repulsion  :  135.87324654 Eh    3697.29901 eV
            if 'Nuclear Repulsion  :' in line:
                if 'nuclear_repulsion_energy' not in self.data['properties'].keys():
                    self.data['properties']['nuclear_repulsion_energy'] = []
                self.data['properties']['nuclear_repulsion_energy'].append(
                    float(line.split()[3])
                )

            # One Electron Energy: -674.26034691 Eh  -18347.55681 eV
            if 'One Electron Energy:' in line:
                if 'scf_one_electron_energy' not in self.data['properties'].keys():
                    self.data['properties']['scf_one_electron_energy'] = []
                self.data['properties']['scf_one_electron_energy'].append(
                    float(line.split()[3])
                )

            # Two Electron Energy:  245.90403408 Eh    6691.38895 eV
            if 'Two Electron Energy:' in line:
                if 'scf_two_electron_energy' not in self.data['properties'].keys():
                    self.data['properties']['scf_two_electron_energy'] = []
                self.data['properties']['scf_two_electron_energy'].append(
                    float(line.split()[3])
                )

            # E(XC)              :      -26.170411238000 Eh 
            if 'E(XC) ' in line:
                if 'scf_xc_energy' not in self.data['properties'].keys():
                    self.data['properties']['scf_xc_energy'] = []
                self.data['properties']['scf_xc_energy'].append(
                    float(line.split()[2])
                )

            line = next(outfile)
    
    def _extract_mp_energies(self, outfile, line):
        """Moller-Plesset calculation properties.

        This is called directly after the ``'ORCA  MP2 '`` trigger, and 
        will terminate once the ``'ORCA property calculations'`` trigger is reached.

        Instead of returning the energies themselves, we handle the creation and
        modification of ``mp_info`` here so any missing information is not an
        issue.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        while '-     ORCA property calculations      *' != line.strip():
            #  MP2 CORRELATION ENERGY   :     -3.132364939 Eh
            if 'MP2 CORRELATION ENERGY' in line:
                if 'mp2_correlation_energy' not in self.data['properties'].keys():
                    self.data['properties']['mp2_correlation_energy'] = []
                if 'RI-MP2' in line:
                    index = 3
                else:
                    index = 4
                self.data['properties']['mp2_correlation_energy'].append(
                    float(line.split()[index])
                )
                break
            
            line = next(outfile)

    def _extract_scf_info(self, outfile, line):
        """Other scf information.

        This will be placed under the ``'keyword'`` JSON property.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.

        Returns
        -------
        :obj:`dict`
            Available SCF information that could contain the following keys.

            ``'rij_approximation'``
                The resolution of identity (RI) approximation for the Coulomb
                (J) term.
            ``'cosx_approximation'``
                The chain-of-spheres integration approximation to the exchange
                term (COSX).
        """
        while 'Total time needed     ' not in line.strip():
            #  RI-approximation to the Coulomb term is turned on
            if 'RI-approximation to the Coulomb term is turned on' in line:
                self.data['keywords']['rij_approximation'] = True

            #    RIJ-COSX (HFX calculated with COS-X)).... on
            if 'RIJ-COSX (HFX calculated with COS-X)' in line:
                self.data['keywords']['cosx_approximation'] = True
            
            if 'RI-JK (J+K treated both via RI)' in line:
                self.data['keywords']['rik_approximation'] = True
            
            line = next(outfile)
    
    def _extract_mulliken_charges(self, outfile, line):
        """Mulliken atomic charges in same order as atomic coordinates.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        line = next(outfile)
        line = next(outfile)

        # Creates initial mulliken_charges property.
        if 'mulliken_charges' not in self.data['properties'].keys():
            self.data['properties']['mulliken_charges'] = []
        
        # Appends Mulliken charges to a new item for every structure.
        self.data['properties']['mulliken_charges'].append([])
        while 'Sum of atomic charges' not in line:
            line_split = line.split(':')
            self.data['properties']['mulliken_charges'][-1].append(
                float(line_split[-1])
            )
            line = next(outfile)
    
    def _extract_loewdin_charges(self, outfile, line):
        """Loewdin atomic charges in same order as atomic coordinates.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        line = next(outfile)
        line = next(outfile)

        # Creates initial loewdin_charges property.
        if 'loewdin_charges' not in self.data['properties'].keys():
            self.data['properties']['loewdin_charges'] = []
        
        # Appends Loewdin charges to a new item for every structure.
        self.data['properties']['loewdin_charges'].append([])
        while '' != line.strip():
            line_split = line.split(':')
            self.data['properties']['loewdin_charges'][-1].append(
                float(line_split[-1])
            )
            line = next(outfile)
    
    def _extract_dipole(self, outfile, line):
        """The X, Y, and Z dipole components.

        Final QCJSON specifies the method of the dipole moment (e.g.,
        ``'scf_dipole_moment'``, ``'mp2_dipole_moment'``). For now, we just
        store it as ``'dipole_moment'``.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        if 'dipole_moment' not in self.data['properties'].keys():
            self.data['properties']['dipole_moment'] = []
        
        while 'Total Dipole Moment    :' not in line:
            line = next(outfile)
        line_split = line.split()
        dipole = [float(line_split[4]), float(line_split[5]), float(line_split[6])]
        self.data['properties']['dipole_moment'].append(dipole)
    
    def _add_geo_conv(self, info_label, line):
        """Parse and add geometric convergence info to data.

        Parameters
        ----------
        info_label : :obj:`str`
            Label for geometric convergence criteria.
        line : :obj:`str`
            Line from output file to extract information from.
        """
        split_line = line.split()
        value = float(split_line[2])
        target = float(split_line[3])
        if f'geo_{info_label}_target' not in self.data['keywords'].keys():
            self.data['keywords'][f'geo_{info_label}_target'] = target
        try:
            self.data['keywords'][f'geo_{info_label}_value'].append(value)
        except KeyError:
            self.data['keywords'][f'geo_{info_label}_value'] = [value]
    
    def _extract_geo_conv(self, outfile, line):
        """Extract geometric convergence values and tolerance.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        while 'Max(Dihed)' not in line and 'Max(Improp)' not in line:
            if 'Energy change' in line:
                self._add_geo_conv('energy_change', line)
            elif 'RMS gradient' in line:
                self._add_geo_conv('rms_gradient', line)
            elif 'MAX gradient' in line:
                self._add_geo_conv('max_gradient', line)
            elif 'RMS step' in line:
                self._add_geo_conv('rms_step', line)
            elif 'MAX step' in line:
                self._add_geo_conv('max_step', line)
            line = next(outfile)
    
    def _extract_frequencies(self, outfile, line):
        """Extract vibrational frequencies, omegas. Includes 0.00 freqencies.

        Based on https://github.com/MolSSI/QCSchema/pull/50#issuecomment-499155251.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        # Skips the following lines.
        line = next(outfile)  # -----------------------
        line = next(outfile)  # 
        line = next(outfile)  # Scaling factor for frequencies =  1.000000000
        line = next(outfile)  # 
        line = next(outfile)  #    0:         0.00 cm**-1

        # Sets up data.
        if 'omega' not in self.data['properties'].keys():
            self.data['properties']['omega'] = []

        vibfreqs = []
        while line.strip() != '':
            freq = float(line.split()[1])
            vibfreqs.append(freq)
            
            line = next(outfile)

        self.data['properties']['omega'].append(vibfreqs)

    def _extract_normal_modes(self, outfile, line):
        """Extract normalized, mass-weighted normal modes, q.

        Based on https://github.com/MolSSI/QCSchema/pull/50#issuecomment-499155251.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        for _ in range(0, 8):
            line = next(outfile)

        # Sets up data.
        if 'q' not in self.data['properties'].keys():
            self.data['properties']['q'] = []

        q = []
        mode = 0
        while line.strip() != '-----------':
            if '.' not in line:
                # Skips lines that print column numbers. For example,
                #   0          1          2          3          4          5 
                # Also restarts mode index.
                mode = 0
            else:
                disp = [float(i) for i in line.split()[1:]]
                try:
                    q[mode].extend(disp)
                except IndexError:
                    q.append(disp)
                mode += 1
            line = next(outfile)
        self.data['properties']['q'].append(q)
    
    def _extract_temp(self, outfile, line):
        """Extracts temperature used for thermochemistry.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        # Sets up data.
        if 'temperature' not in self.data['properties'].keys():
            self.data['properties']['temperature'] = []
        temp = float(line.split()[2])
        self.data['properties']['temperature'].append(temp)
        # Needs to move off this line so extract can continue.
        line = next(outfile)
    
    def _extract_zpve(self, outfile, line):
        """Extract zero-point vibrational energy correction.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        # Sets up data.
        if 'zero_point_vibrational_correction' not in self.data['properties'].keys():
            self.data['properties']['zero_point_vibrational_correction'] = []
        zpve = float(line.split()[4])
        self.data['properties']['zero_point_vibrational_correction'].append(zpve)
        # Needs to move off this line so extract can continue.
        line = next(outfile)
    
    def _extract_thermal_corrections(self, outfile, line):
        """Extracts thermal vibrational, rotational, and translational
        corrections.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        # Sets up data.
        if 'thermal_energy_corrections' not in self.data['properties'].keys():
            self.data['properties']['thermal_energy_corrections'] = []
        thermal = float(line.split()[3])
        self.data['properties']['thermal_energy_corrections'].append(thermal)
        # Needs to move off this line so extract can continue.
        line = next(outfile)
    
    def _extract_enthalpic_corr(self, outfile, line):
        """Extract enthalpic corrections; the difference between H_298 and
        E_298.

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        # Sets up data.
        if 'enthalpic_corrections' not in self.data['properties'].keys():
            self.data['properties']['enthalpic_corrections'] = []
        enthalpy = float(line.split()[4])
        self.data['properties']['enthalpic_corrections'].append(enthalpy)
        # Needs to move off this line so extract can continue.
        line = next(outfile)
    
    def _extract_entropic_corr(self, outfile, line):
        """Extract translational, rotational, and vibrational entropic
        corrections to enthalpy for Gibbs free energy (i.e., T*S_298).

        Parameters
        ----------
        outfile : :obj:`io.TextIOWrapper`
            Buffered text stream of the output file.
        line : :obj:`str`
            Parsed line from ``outfile``.
        """
        # Sets up data.
        if 'entropic_corrections' not in self.data['properties'].keys():
            self.data['properties']['entropic_corrections'] = []
        entropy = float(line.split()[4])
        self.data['properties']['entropic_corrections'].append(entropy)
        # Needs to move off this line so extract can continue.
        line = next(outfile)

