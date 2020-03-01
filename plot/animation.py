# This file is part of CELADRO, Copyright (C) 2016-20, Romain Mueller
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
import matplotlib.pyplot as plt
import matplotlib.animation as ani


######################################################################
# animation


def animate(oa, fn, rng=[], inter=200, show=True):
    """
    Show a frame-by-frame animation.

    Args:
        oa -- the output archive
        fn -- the plot function (argument: frame, plot engine)
        rng -- range of the frames to be ploted
        interval -- time between frames (ms)
    """
    # set range
    if len(rng) == 0:
        rng = [1, oa._nframes+1]

    # create the figure
    fig = plt.figure()

    # the local animation function
    def animate_fn(i):
        # we want a fresh figure everytime
        fig.clf()
        # load the frame
        frame = oa.read_frame(i)
        # call the global function
        fn(frame, fig)

    anim = ani.FuncAnimation(fig, animate_fn,
                             frames=np.arange(rng[0], rng[1]),
                             interval=inter, blit=False)

    if show:
        return plt.show()
    else:
        return anim


def save(an, fname, fps, tt='ffmpeg', bitrate=-1):
    writer = ani.writers[tt](fps=fps, bitrate=bitrate)
    an.save(fname, writer=writer)
