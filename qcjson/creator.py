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

from qcjson import __version__ as qcjson_version
from qcjson.utils import get_files
from qcjson.utils import select_files
from qcjson.utils import cclib_version_check

from qcjson.jsons import orcaJSON
from qcjson.jsons import xtbJSON

def identify_package(outfile_path):
    """Identifies computational chemistry package.

    Only supported packaged should be included in `triggers`.

    Parameters
    ----------
    outfile_path : :obj:`str`
        Path to output file.
    
    Returns
    -------
    :obj:`chemicalJSON`
        One of the supported ChemicalJSON classes.
    """
    with open(outfile_path, 'r') as f:
        for line in f:
            for parser, phrases, do_break in triggers:
                if all([line.lower().find(p.lower()) >= 0 for p in phrases]):
                    filetype = parser
                    if do_break:
                        return filetype

### Runtime Functions ###

# Triggers to identify output files.
triggers = [
    (orcaJSON, ["O   R   C   A"], True),
    (xtbJSON, ["x T B"], True)
]

def qcjson_creator(output_file, save_dir, debug, prettify):
    """Creates a single QCJSON file.

    This is called by the qcjson-creator.py script.

    Parameters
    ----------
    output_file : :obj:`str`
        Path to output file.
    save_dir : :obj:`str`
        Directory to save the QCJSON to.
    debug : :obj:`bool`
        Whether or not to raise errors during the QCJSON process instead of
        just skipping over the file.
    prettify : :obj:`bool`
        Indent each JSON property.
    
    Returns
    -------
    :obj:`list`
        All file paths that encountered errors if ``debug`` is ``True``.
    """
    geom_splits = [xtbJSON]

    json_package = identify_package(output_file)
    
    if json_package in geom_splits:
        raise TypeError('This code needs to be used as geom_split')
    
    out_json = json_package(output_file)
    json_dict = out_json.get_json(debug=debug)
    if out_json.path not in out_json.error_files:
        if save_dir == './':
            abs_path = os.path.dirname(out_json.path)
        else:
            abs_path = save_dir
        
        out_json.write(
            out_json.name, json_dict, abs_path, prettify=prettify
        )
    
    return out_json.error_files

def qcjson_creator_split(output_file, geom_file, save_dir, debug, prettify):
    """Creates a single QCJSON file from an output file and xyz file.

    This is called by the qcjson-creator.py script.

    Parameters
    ----------
    output_file : :obj:`str`
        Path to output file.
    geom_file : :obj:`str`
        Path to file containing xyz geometries.
    save_dir : :obj:`str`
        Directory to save the QCJSON to.
    debug : :obj:`bool`
        Whether or not to raise errors during the QCJSON process instead of
        just skipping over the file.
    prettify : :obj:`bool`
        Indent each JSON property.
    
    Returns
    -------
    :obj:`list`
        All file paths that encountered errors if ``debug`` is ``True``.
    """
    # pylint: disable=too-many-function-args
    json_package = identify_package(output_file)
    out_json = json_package(output_file, geom_file)
    json_dict = out_json.get_json(debug=debug)
    if out_json.path not in out_json.error_files:
        if save_dir == './':
            abs_path = os.path.dirname(out_json.path)
        else:
            abs_path = save_dir
        
        out_json.write(
            out_json.name, json_dict, abs_path, prettify=prettify
        )

    return out_json.error_files