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
import argparse

from qcjson import __version__ as qcjson_version
from qcjson.creator import qcjson_creator
from qcjson.creator import qcjson_creator_split
from qcjson.utils import cclib_version_check
from qcjson.utils import get_files
from qcjson.utils import select_files

def main():
    
    parser = argparse.ArgumentParser(
        description='Creates QCJSONs from output files using QCSchemas. Files must '
                    'have "out" somewhere in the file name or as an extension.'
    )
    parser.add_argument(
        'outputs', metavar='outputs', type=str, nargs='?', default='.',
        help='Path to directory or specific computational chemistry output file.'
    )
    parser.add_argument(
        '--geom_split', metavar='geom_split', type=str, nargs='?', default='',
        help='XYZ file containing all structures from the calculation.'
    )
    parser.add_argument(
        '--save_dir', metavar='save_dir', type=str, nargs='?', default='.',
        help='Path to save JSON files.'
    )
    parser.add_argument(
        '-r', '--recursive', action='store_true',
        help='Recursively create JSONs.'
    )
    parser.add_argument(
        '-o', '--overwrite', action='store_true', help='Overwrite JSON files'
    )
    parser.add_argument(
        '-p', '--prettify', action='store_true',
        help='Prettify JSON files with indentation.'
    )
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help='Will not continue if an error is encountered'
    )
    parser.add_argument(
        '--exclude', nargs='+', default=[],
        help='Ignore paths that contain at least one of these words.'
    )
    parser.add_argument(
        '--include', nargs='+', default=[],
        help='Only include files that contain all of these words.'
    )

    args = parser.parse_args()

    print(f'QCJSON creator v{qcjson_version}')
    print('Written by Alex M. Maldonado (@aalexmmaldonado)')
    print('Energies and distances are Hartrees and Angstroms\n')

    cclib_version_check()

    save_dir = args.save_dir
    if save_dir[-1] != '/':
        save_dir += '/'
    outputs = args.outputs
    geomfile = args.geom_split

    if geomfile != '':
        if args.recursive:
            print('Split geometry output files are not supported recursively.')
        
        print(f'Making QCJSON for {outputs}')
        error_files = qcjson_creator_split(
            outputs, geomfile, save_dir, args.debug, args.prettify
        )

        print(f'\n{len(error_files)} file(s) encountered errors and not written')
        for i in error_files:
            i_name = os.path.basename(i)
            print(f'\u001b[31;1m    {i_name}\u001b[0m')
    else:
        all_error_files = []

        # A file was provided for the outputs.
        if os.path.isfile(outputs):
            print(f'Making QCJSON for {outputs}')
            qcjson_creator(outputs, save_dir, args.debug, args.prettify)
        
        # A directory was provided for the outputs.
        elif os.path.isdir(outputs):

            if outputs[-1] != '/':
                outputs += '/'
            
            if args.recursive:
                print(f'Looking for output files in {outputs}, recursively')
            else:
                print(f'Looking for output files in {outputs}')
            all_outfiles = get_files(outputs, 'out', recursive=args.recursive)
            
            print(f'Found {len(all_outfiles)} output files\n')

            all_outfiles = select_files(all_outfiles, exclude=args.exclude, include=args.include)

            for outfile in all_outfiles:
                file_name = '.'.join(os.path.basename(outfile).split('.')[:-1])
                if save_dir == './':
                    abs_path = os.path.dirname(os.path.abspath(outfile))
                    if not args.overwrite \
                    and os.path.exists(f'{abs_path}/{file_name}.json'):
                        print(
                            f'\n\u001b[36;1m{file_name}.json already exists.\u001b[0m'
                        )
                        continue
                else:
                    if not args.overwrite \
                    and os.path.exists(f'{save_dir}/{file_name}.json'):
                        print(
                            f'\n\u001b[36;1m{file_name}.json already exists.\u001b[0m'
                        )
                        continue
                print(f'Making QCJSON for {file_name}')
                error_files = qcjson_creator(
                    outfile, save_dir, args.debug, args.prettify
                )
                all_error_files.extend(error_files)
                
        else:
            raise ValueError(f'{outputs} is an unsupported type.')
        
        print(f'\n{len(all_error_files)} file(s) encountered errors and not written')
        for i in all_error_files:
            i_name = os.path.basename(i)
            print(f'\u001b[31;1m    {i_name}\u001b[0m')


if __name__ == "__main__":
    main()
