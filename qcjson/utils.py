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
import json
from packaging import version

### Hard coded information ###

_element_to_z = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9,
    'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16,
    'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23,
    'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
    'Ga': 31,'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
    'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44,
    'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51,
    'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58,
    'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65,
    'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72,
    'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
    'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86,
    'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93,
    'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100,
    'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106,
    'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112,
    'Uuq': 114, 'Uuh': 116,
}
_z_to_element = {v: k for k, v in _element_to_z.items()}

# Program-specific keywords categorized into scf, moller-plesset, and coupled
# cluster methods as specified by QCSchema. All methods are in lowercase to 
# assist identification.
methods = {
    'orca': {
        'scf': [
            'rhf', 'rks', 'uhf', 'uks', 'rohf', 'roks', 'hf',
            'hfs', 'lda', 'lsd', 'vwn', 'vwn5', 'vwn3', 'pwlda', 'bp86', 'bp',
            'blyp', 'olyp', 'glyp', 'xlyp', 'pw91', 'mpwpw', 'mpwlyp', 'mpwlyp',
            'pbe', 'rpbe', 'revpbe', 'pwp', 'b1lyp', 'b3lyp', 'b3lyp/g',
            'o3lyp', 'x3lyp', 'b1p', 'b3p', 'b3pw', 'pw1pw', 'mpw1pw',
            'mpw1lyp', 'pbe0', 'pw6b95', 'bhandhylpl', 'tpss', 'tpssh', 'tpss0',
            'm06l', 'mo6', 'm062x', 'b97m-v', 'b97m-d3bj', 'scanfunc',
            'wb97', 'wb97x', 'wb97x-d3', 'wb97x-v', 'wb97x-d3bj', 'cam-b3lyp',
            'lc-blyp', 'b2plyp', 'b2plyp-d', 'b2plyp-d3', 'mpw2plpy',
            'mpw2plyp-d', 'b2gp-plyp', 'b2k-plyp', 'b2t-plyp', 'pwpb95',
            'dsd-blyp', 'dsd-pbep86', 'dsd-pbep95', 'wb2plyp', 'wb2gp-plyp'
        ],
        'moller-plesset': [
            'mp2', 'ri-mp2', 'scs-mp2', 'ri-scs-mp2', 'oo-ri-mp2',
            'oo-ri-scs-mp2', 'mp2-f12', 'mp2-f12', 'mp2-f12-ri', 'mp2-f12d-ri',
            'mp3', 'scs-mp3', 'dlpno-mp2', 'dlpno-scs-mp2', 'dlpno-mp2-f12',
            'dlpno-mp2-f12/d'
        ],
        'coupled cluster': [
            'ccsd', 'ccsd(t)', 'ccsd-f12', 'ccsd(t)-f12', 'ccsd-f12/ri',
            'ccsd(t)-f12d/ri', 'lpno-ccsd', 'dlpno-ccsd', 'dlpno-ccsd(t)',
            'dlpno-ccsd(t1)', 'dlpno-ccsd-f12', 'dlpno-ccsd-f12/d'
        ]
    }
}

# Program-specific, basis-set keywords in lowercase (to help with
# identification).
basis_sets = {
    'orca': [
        '3-21g', 'sto-3g', '3-21gsp', '4-22gsp', '6-31g', 'm6-31g', '6-311g',
        '6-31g*', '6-31g(d)', '6-31g**', '6-31g(d,p)', '6-31g(2d)',
        '6-31g(2df)', '6-31g(2d,p)', '6-31g(2d,2p)', '6-31g(2df,2dp)',
        '6-311g*', '6-311g(d)', '6-311g**', '6-311g(d,p)', '6-311g(2d)',
        '6-311g(2df)', '6-311g(2d,p)', '6-311g(2d,2p)', '6-311g(2df,2dp)',
        '6-311g(3df)','6-311g(3df,3pd)', '6-31+g', '6-31++g(d,p)',
        'def2-svp', 'def2-sv(p)', 'def2-tzvp', 'def2-tzvp(-f)', 'def2-tzvpp',
        'def2-qzvpp', 'sv', 'sv(p)', 'svp', 'tzv', 'tzv(p)', 'tzvp', 'tzvpp',
        'qzvp', 'qzvpp', 'ma-def2-svp', 'ma-def2-sv(p)', 'ma-def2-tavp',
        'ma-def2-tzvp(-f)', 'ma-def2-tzvpp', 'ma-def2-qzvpp', 'def2-svpd',
        'def2-tzvpd', 'def2-tzvppd', 'def2-qzvpd', 'def2-qzvppd',
        'dkh-def2-svp', 'zora-def2-svp', 'dkh-def2-sv(p)', 'zora-def2-sv(p)',
        'dkh-def2-tzvp', 'zora-def2-tzvp', 'dhk-def2-tzvp(-f)',
        'zora-def2-tzvp(-f)', 'dhk-def2-tzvpp', 'zora-def2-tzvpp',
        'dkh-def2-qzvpp', 'zora-def2-qzvpp', 'ma-dkh-def2-svp',
        'ma-zora-def2-svp', 'ma-dkh-def2-sv(p)', 'ma-zora-def2-sv(p)',
        'ma-dkh-def2-tzvp', 'ma-zora-def2-tzvp', 'ma-dhk-def2-tzvp(-f)',
        'ma-zora-def2-tzvp(-f)', 'ma-dhk-def2-tzvpp', 'ma-zora-def2-tzvpp',
        'ma-dkh-def2-qzvpp', 'ma-zora-def2-qzvpp', 'dkh-svp', 'zora-svp',
        'dkh-sv(p)', 'zora-sv(p)', 'dkh-tzvp', 'zora-tzvp', 'dhk-tzvp(-f)',
        'zora-tzvp(-f)', 'dhk-tzvpp', 'zora-tzvpp', 'dkh-qzvpp', 'zora-qzvpp',
        'ma-dkh-svp', 'ma-zora-svp', 'ma-dkh-sv(p)', 'ma-zora-sv(p)',
        'ma-dkh-tzvp', 'ma-zora-tzvp', 'ma-dhk-tzvp(-f)', 'ma-zora-tzvp(-f)',
        'ma-dhk-tzvpp', 'ma-zora-tzvpp', 'ma-dkh-qzvpp', 'ma-zora-qzvpp',
        'sarc-dkh-tzvp', 'sarc-dkh-tzvpp', 'sarc-zora-tzvp', 'sarc-zora-tzvpp',
        'sarc-dkh-svp', 'sarc-zora-svp', 'sarc2-dkh-qzv', 'sarc2-dkh-qzvp',
        'sarc2-zora-qzvp', 'sarc-zora-qzvp', 'pc-0' ,'pc-1', 'pc-2', 'pc-3',
        'pc-4', 'aug-pc-0' ,'aug-pc-1', 'aug-pc-2', 'aug-pc-3', 'aug-pc-4',
        'pcseg-0' ,'pcseg-1', 'pcseg-2', 'pcseg-3', 'pcseg-4', 'aug-pcseg-0',
        'aug-pcseg-1', 'aug-pcseg-2', 'aug-pcseg-3', 'aug-pcseg-4', 'pcsseg-0',
        'pcsseg-1', 'pcsseg-2', 'pcsseg-3', 'pcsseg-4', 'aug-pcsseg-0',
        'aug-pcsseg-1', 'aug-pcsseg-2', 'aug-pcsseg-3', 'aug-pcsseg-4', 'pcj-0',
        'pcj-1', 'pcj-2', 'pcj-3', 'pcj-4', 'aug-pcj-0', 'aug-pcj-1',
        'aug-pcj-2', 'aug-pcj-3', 'aug-pcj-4', 'sapporo-dzp-2012',
        'sapporo-tzp-2012', 'sapporo-qzp-2012', 'sapporo-dkh3-dzp-2012',
        'sapporo-dkh3-tzp-2012', 'sapporo-dkh3-qzp-2012', 'cc-pvdz', 'cc-pvtz',
        'cc-pvqz', 'cc-pv5z', 'cc-pv6z', 'aug-cc-pvdz', 'aug-cc-pvtz',
        'aug-cc-pvqz', 'aug-cc-pv5z', 'aug-cc-pv6z', 'cc-pcvdz', 'cc-pcvtz',
        'cc-pcvqz' 'cc-pcv5z', 'cc-pcv6z', 'aug-cc-pcvdz', 'aug-cc-pcvtz',
        'aug-cc-pcvqz' 'aug-cc-pcv5z', 'aug-cc-pcv6z', 'cc-pwcvdz', 'cc-pwcvtz',
        'cc-pwcvqz' 'cc-pwcv5z', 'aug-cc-pwcvdz', 'aug-cc-pwcvtz',
        'aug-cc-pwcvqz' 'aug-cc-pwcv5z', 'aug-cc-pwcvd(+d)z',
        'aug-cc-pwcvt(+d)z', 'aug-cc-pwcvq(+d)z' 'aug-cc-pwcv5(+d)z',
        'ano-pvdz', 'ano-pvtz', 'ano-pvqz', 'ano-pv5z', 'saug-ano-pvdz',
        'saug-ano-pvtz', 'saug-ano-pvqz', 'aug-ano-pvdz', 'aug-ano-pvtz',
        'aug-ano-pvqz', 'ano-rcc-full', 'ano-rcc-dzp', 'ano-rcc-tzp',
        'ano-rcc-qzp', 'd95', 'd95p', 'mini', 'minis', 'midi', 'minix',
        'wachters+f', 'partridge-1', 'partridge-2', 'partridge-3',
        'partridge-4', 'lanl2dz', 'lanl2tz', 'lanl2tz(f)', 'lanl08', 'lanl08(f)',
        'epr-ii', 'epr-iii', 'iglo-ii', 'iglo-iii', 'aug-cc-pvtz-j'
    ]
}

### Utility functions ###

def atoms_by_element(atom_list):
    """Converts a list of atoms identified by their atomic number to their
    elemental symbol in the same order.
    
    Parameters
    ----------
    atom_list : :obj:`list` [:obj:`int`]
        Atomic numbers of atoms within a structure.
    
    Returns
    -------
    :obj:`list` [:obj:`str`]
        Element symbols of atoms within a structure.
    """
    return [_z_to_element[i] for i in atom_list]

def get_files(path, expression, recursive=True):
    """Returns paths to all files in a given directory that matches a provided
    expression in the file name. Commonly used to find all files of a certain
    type, e.g. output or xyz files.
    
    Parameters
    ----------
    path : :obj:`str`
        Specifies the directory to search.
    expression : :obj:`str`
        Expression to be tested against all file names in 'path'.
    recursive :obj:`bool`, optional
        Recursively find all files in all subdirectories.
    
    Returns
    -------
    :obj:`list` [:obj:`str`]
        All absolute paths to files matching the provided expression.
    """
    if path[-1] != '/':
        path += '/'
    if recursive:
        all_files = []
        for (dirpath, _, filenames) in os.walk(path):
            index = 0
            while index < len(filenames):
                if dirpath[-1] != '/':
                    dirpath += '/'
                filenames[index] = dirpath + filenames[index]
                index += 1
            all_files.extend(filenames)
        files = []
        for f in all_files:
            if expression in f:
                files.append(f)
    else:
        files = []
        for f in os.listdir(path):
            filename = os.path.basename(f)
            if expression in filename:
                files.append(path + f)
    return files

def convert_forces(
    forces, e_units_calc, r_units_calc, e_units, r_units
):
    """Converts forces (or gradients) to specified units.

    Parameters
    ----------
    forces : :obj:`numpy.ndarray`
        An array with units of energy and distance matching `e_units_calc`
        and `r_units_calc`.
    e_units_calc : :obj:`str`
        Specifies package-specific energy units used in calculation. Available
        units are ``'eV'``, ``'hartree'``, ``'kcal/mol'``, and ``'kJ/mol'``.
    r_units_calc : :obj:`str`
        Specifies package-specific distance units used in calculation. Available
        units are ``'Angstrom'`` and ``'bohr'``.
    e_units : :obj:`str`
        Desired units of energy. Available units are ``'eV'``, ``'hartree'``,
        ``'kcal/mol'``, and ``'kJ/mol'``.
    r_units : obj:`str`
        Desired units of distance. Available units are ``'Angstrom'`` and
        ``'bohr'``.
    
    Returns
    -------
    :obj:`numpy.ndarray`
        Forces converted into the desired units.
    """
    #'ORCA': {'e_unit': 'hartree', 'r_unit': 'bohr'}
    if e_units not in ['eV', 'hartree', 'kcal/mol', 'kJ/mol']:
        raise ValueError(f'{e_units} is not an available energy unit.')
    if r_units not in ['Angstrom', 'bohr']:
        raise ValueError(f'{r_units} is not an available distance unit.')
    forces_conv = forces
    if e_units_calc != e_units:
        forces_conv = cclib.parser.utils.convertor(
            forces_conv, e_units_calc, e_units
        )
    if r_units_calc != r_units:
        forces_conv = cclib.parser.utils.convertor(
            forces_conv, r_units, r_units_calc
        )
    return forces_conv

def read_json(json_path):
    """Read JSON file.
    
    Parameters
    ----------
    json_path : :obj:`str`
        Path to json file.
    
    Returns
    -------
    :obj:`dict`
        Contents of JSON file.
    """
    with open(json_path, 'r') as reader:
        json_dict = json.loads(reader.readline())
    
    return json_dict

def select_files(all_files, exclude=[], include=[]):
    """Removes all undesired files from list of paths.

    Parameters
    ----------
    all_files : :obj:`list` [:obj:`str`]
        Collection of file paths.
    exclude : :obj:`list` [:obj:`str`]
        Ignore paths that contain at least one of these strings.
    include : :obj:`list` [:obj:`str`]
        Only include files that contain all of these strings.
    
    Returns
    -------
    :obj:`list`
        Paths that are meet the inclusion and exclusion critera.
    """
    # Removes files that match any of the words in the remove list.
    start_number = len(all_files)
    if len(exclude) > 0:
        for trigger in exclude:
            print(f'Selecting files including: {trigger}')
            all_files = [i for i in all_files if trigger not in i]
        end_number = len(all_files)
    if len(include) > 0:
        for trigger in include:
            print(f'Selecting files not including: {trigger}')
            all_files = [i for i in all_files if trigger in i]
        end_number = len(all_files)
    if 'end_number' in locals():
        print(
            f'Removed {start_number-end_number} file(s); {end_number} remain'
        )
    return all_files

### Runtime Functions ###

def cclib_version_check():
    """Ensures cclib version is at least 1.7.

    Dispersion energies were not systematically included or parsed until version
    1.6.4. However, the __version__ property was incorrect until version 1.7.
    """
    cclib_version = cclib.__version__
    if version.parse(cclib_version) < version.parse('1.7'):
        raise ValueError(
            f'cclib version is {cclib_version}; need 1.7 or higher.'
        )
