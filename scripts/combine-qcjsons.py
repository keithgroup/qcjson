#!/usr/bin/env python3

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
from functools import reduce
import argparse

from qcjson import __version__ as qcjson_version
from qcjson.jsons.qcjson import QCJSON
from qcjson.utils import read_json

def dict_cleaner(data):
    new_data = {}
    for k, v in data.items():
        if isinstance(v, dict):
            v = dict_cleaner(v)
        if not v in (u'', None, {}):
            new_data[k] = v
    return new_data

def combine_jsons(json_dir, nest_jsons):
    """Writes a combined JSON file containing all files.

    Parameters
    ----------
    json_dir : :obj:`str`
        Path to directory where chemical JSON files will be saved.
    name : :obj:`str`
        File name for cumulative chemical JSON file.
    """
    all_json_dict = {}

    start = json_dir.rfind(os.sep) + 1
    # Walks through every dir inside save_dir
    for path, _, files in os.walk(json_dir):

        files = [i for i in files if '.json' in i]  # List of json files in current search dir.
        folders = path[start:].split(os.sep)  # List of the folders from the save dir to current search dir.
        if path[-1] != '/':
            path += '/'
        try:
            json_files = {}
            for json_file in files:
                json_path = path + json_file
                json_dict = read_json(json_path)
                if isinstance(json_dict, list):
                    json_name = json_dict[0]['name']
                else:
                    json_name = json_dict['name']
                json_files[json_name] = json_dict

            if len(json_files) == 0:
                # Creates empty dicts for each dir.
                parent = reduce(dict.get, folders[:-1], all_json_dict)
                parent[folders[-1]] = json_files
            elif len(json_files) == 1:
                if not nest_jsons:
                    # A directory only has one JSON file. Most likely a single
                    # calculation folder and should not be a nested dict with only
                    # one value. Going to replace this dir dict with this single 
                    # JSON file.
                    parent = reduce(dict.get, folders[:-1], all_json_dict)
                    key = list(json_files.keys())[0]
                    parent[key] = json_files[key]
                else:
                    # Will have nested dictionaries if there is one json per dir.
                    parent = reduce(dict.get, folders[:-1], all_json_dict)
                    parent[folders[-1]] = json_files
            else:
                parent = reduce(dict.get, folders[:-1], all_json_dict)
                parent[folders[-1]] = json_files
        except KeyError:
            # Will not include JSON files without 'name' key.
            # Useful for using --overwrite when cumulative file is present.
            pass
    
    # Checks and removes empty dictionaries.
    all_json_dict = dict_cleaner(all_json_dict)

    return all_json_dict

def main():
    
    parser = argparse.ArgumentParser(
        description='Combines all JSON files into a single JSON organized by '
                    'their directories. Files must end in "json".'
    )
    parser.add_argument(
        'json_dir', metavar='json_dir', type=str, nargs='?', default='.',
        help='Path to directory containing JSON files.'
    )
    parser.add_argument(
        '--save_dir', metavar='save_dir', type=str, nargs='?', default='.',
        help='Path to save the combined JSON file.'
    )
    parser.add_argument(
        '-n', '--name', metavar='name', type=str, nargs='?', default='combined',
        help='Name of the combined JSON file.'
    )
    parser.add_argument(
        '--nest_jsons', action='store_true', help='Dirs with one JSON will be nested.'
    )
    parser.add_argument(
        '-r', '--recursive', action='store_true',
        help='Recursively create JSONs.'
    )
    parser.add_argument(
        '-o', '--overwrite', action='store_true', help='Overwrite JSON files.'
    )
    parser.add_argument(
        '-p', '--prettify', action='store_true',
        help='Prettify JSON files with indentation.'
    )

    args = parser.parse_args()

    print(f'Combine QCJSONs v{qcjson_version}')
    print('Written by Alex M. Maldonado (@aalexmmaldonado)')
    print('Energies and distances are Hartrees and Angstroms\n')

    json_dir = os.path.abspath(args.json_dir)
    save_dir = os.path.abspath(args.save_dir)
    if json_dir[-1] != '/':
        json_dir += '/'
    if save_dir[-1] != '/':
        save_dir += '/'

    if os.path.exists(save_dir + args.name + '.json') and not args.overwrite:
        print('ERROR: Combined JSON file already exists and overwrite is False.')
        print('Will not continue.')
        exit()
    
    # Gets all JSON file paths.
    looking_string = f'Looking for JSON files in {json_dir}'
    if args.recursive:
        looking_string += ', recursively'
    print(looking_string)
    
    json_dict = combine_jsons(json_dir, args.nest_jsons)

    qcjson = QCJSON()
    qcjson.write(args.name, json_dict, save_dir, prettify=args.prettify)







if __name__ == "__main__":
    main()
