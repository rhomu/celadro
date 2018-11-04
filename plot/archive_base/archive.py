# This file is part of CELADRO, Copyright (C) 2016-17, Romain Mueller
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
from zipfile import ZipFile
import json
import numpy as np


class archive:
    """
    Import archive.

    Automatically import the parameters and allows to extract the individual
    frame files as well.
    """
    def __init__(self, path):
        """Reads parameters from archive."""
        # get extension and file name
        self._path = path
        self._name, ext = os.path.splitext(path)
        self._ext = '.json'

        # check if we use a compressed archive or a directory
        if ext == "":
            self._compress_full = False
        elif ext == '.zip':
            self._compress_full = True
        else:
            raise ValueError('Archive type ' + path + ' not recognized')

        # check for single file compression
        if os.path.isfile(os.path.join(self._path, 'parameters.json.zip')):
            self._compress = True
            self._ext = '.json.zip'
        else:
            self._compress = False

        # load the parameters as a dictionary
        dat = self.extract_and_read('parameters')
        self.parameters = {entry: self.get_value(dat[entry]['value'],
                           dat[entry]['type']) for entry in dat}
        # add the variables to the object (i luv python)
        self.__dict__.update(self.parameters)
        # total number of frames
        self._nframes = int((self.nsteps-self.nstart)/self.ninfo)

    def get_value(self, v, t):
        """Convert string to value with correct type handling."""
        if v == "nan" or v == "-nan":
            raise ValueError('Nan found while converting to ' + t)

        if t == 'double' or t == 'float':
            return float(v)
        elif t == 'int' or t == 'unsigned':
            return int(v)
        elif t == 'long' or t == 'unsigned long':
            return int(v)
        elif t == 'bool':
            return bool(v)
        elif t == 'string':
            # quotes are already erased for string types
            return v
        elif t[:5] == 'array':
            return np.array([self.get_value(i, t[6:-1]) for i in v])
        else:
            raise ValueError('Unrecognized type ' + t)

    def extract_and_read(self, fname):
        """Extract json file from archive."""
        # extract
        if self._compress_full:
            with ZipFile(self._path, "r") as f:
                data = f.read(fname+self._ext)
            return json.loads(data)['data']
        elif self._compress:
            with ZipFile(os.path.join(self._path, fname + self._ext)) as f:
                data = f.read(fname+'.json')
            return json.loads(data)['data']
        else:
            # read content of file
            output = open(os.path.join(self._path, fname + self._ext))
            # load json
            return json.load(output)['data']

    def read_frame(self, frame):
        """Read state file from archive.

        Parameters:
        frame -- the frame number to be read (0,1,2...)
        """
        if frame > self._nframes:
            raise ValueError('Frame does not exist.')
        # get the json
        dat = self.extract_and_read('frame' + str(self.nstart+frame*self.ninfo))
        # convert to dict
        dat = {entry: self.get_value(dat[entry]['value'], dat[entry]['type']) for entry in dat}

        # return a dummy class that holds the data
        class frame_holder:
            """
            Dummy frame holder.

            Automatically define all the variables defined in the corresponding
            json file.
            """
            def __init__(self, parameters):
                self.parameters = parameters
        # create holder and forward parameters
        frame = frame_holder(self.parameters)
        frame.__dict__.update(dat)
        return frame

    def __getitem__(self, frame):
        return self.read_frame(frame)

    def read_frames(self):
        """Generates all frames successively"""
        for n in range(self._nframes+1):
            yield self.read_frame(n)


def loadarchive(path):
    return archive(path)
