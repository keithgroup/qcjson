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

# pylint: disable=no-member

class outfileParser:
    """Base class for parsing output files.

    Parameters
    ----------
    outfile_path : :obj:`str`
        Path to output file.
    
    Attributes
    ----------
    outfile_path : :obj:`str`
        Path to output file.
    file_name : :obj:`str`
        The name of the file without extension.
    data : :obj:`dict`
        Extracted data from the output file (after being parsed).
    """

    def __init__(self, outfile_path):
        self.outfile_path = outfile_path
        self.file_name = '.'.join(outfile_path.split('/')[-1].split('.')[:-1])
        self.data = {
            'keywords': {},
            'properties': {},
            'model': {}
        }
    
    def parse(self):
        """Parses trajectory file and extracts information.
        """
        with open(self.outfile_path, mode='r') as trajfile:
            for trajline in trajfile:
                self.extract(trajfile, trajline)
        self._after_parse()
