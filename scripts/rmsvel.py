#
# This is a simple example file to show the plotting capabilities of the
# program. Uses python2.
#
# Usage:
#
#   python2 plot-cells.py input [output]
#
#  where
#
#    intput -- the input file or directory
#    output -- (optional) if present saves the animation as a video to show to
#              your mom.

import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import numpy as np
from math import sqrt

# import local libs
sys.path.insert(0, "../plot/")
import plot
import archive
import animation

##################################################
# Init

if len(sys.argv)==1:
    print "Please provide an input file."
    exit(1)

# load archive from file
ar = archive.loadarchive(sys.argv[1])

##################################################
# plot simple animation of phases

ar = archive.loadarchive(sys.argv[1])

rms = np.zeros(ar._nframes+1)
tox = np.zeros(ar._nframes+1)
toy = np.zeros(ar._nframes+1)

c = 0
for i in range(0, ar._nframes+1):
    frame = ar.read_frame(i)
    print "{}/{}".format(i, ar._nframes),
    vx, vy = plot.get_velocity_field(frame.phi, frame.velp + frame.velc + frame.velf, size=24)
    rms[c] = sqrt(np.mean(vx**2+vy**2))
    tox[c] = np.mean(vx)
    toy[c] = np.mean(vy)
    print "vx={}, vy={}, rms={}".format(tox[c], toy[c], rms[c])
    c     += 1

plt.plot(np.arange(0, len(rms)), rms)
plt.plot(np.arange(0, len(rms)), tox)
plt.plot(np.arange(0, len(rms)), toy)
plt.show()
