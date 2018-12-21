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

import numpy as np
import archive_base


class archive(archive_base.archive):
    """Simply reshape 2d fields after importing"""

    def __init__(self, path):
        archive_base.archive.__init__(self, path)
        self.parameters['walls'] = np.reshape(self.parameters['walls'],
                                              self.Size)
        self.__dict__.update(self.parameters)

    def read_frame(self, frame):
        frame = archive_base.archive.read_frame(self, frame)

        # array sizes
        lx, ly = self.Size
        px, py = self.patch_size

        phi = []
        for i in range(len(frame.phi)):
            p = np.reshape(frame.phi[i], (px, py))
            # compensate for offset
            p = np.roll(p, frame.offset[i][0], axis=0)
            p = np.roll(p, frame.offset[i][1], axis=1)
            # extend to full size
            p = np.concatenate((p, np.zeros((px, ly-py))), axis=1)
            p = np.concatenate((p, np.zeros((lx-px, ly))), axis=0)
            # put in right postition
            p = np.roll(p, frame.patch_min[i][0], axis=0)
            p = np.roll(p, frame.patch_min[i][1], axis=1)
            # save
            phi.append(p)
        frame.phi = phi

        frame.stress_xx.shape = (lx, ly)
        frame.stress_xy.shape = (lx, ly)
        frame.stress_yy.shape = (lx, ly)

        return frame


def loadarchive(path):
    return archive(path)
