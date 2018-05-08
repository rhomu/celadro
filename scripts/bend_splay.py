#
# This is a simple example file to show the plotting capabilities of the
# program. Uses python2.
#
# Usage:
#
#   python2 bend_splay.py input
#
#  where
#
#    intput -- the input file or directory

import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import numpy as np
from math import sqrt, atan2

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

ll = 50#ar._nframes+1
ar = archive.loadarchive(sys.argv[1])

bend  = np.zeros(ll)
splay = np.zeros(ll)

c = 0
for i in range(0, ll):
    frame = ar.read_frame(i)
    print "{}/{}".format(i, ar._nframes),
    Qxx, Qxy = plot.get_Qtensor(frame.phi, frame.Q00, frame.Q01, size=24)
    S  = np.vectorize(sqrt)(Qxy**2 + Qxx**2)
    nx = np.vectorize(sqrt)((1 + Qxx/S)/2)
    ny = np.sign(Qxy)*np.vectorize(sqrt)((1 - Qxx/S)/2)
    nn = nx*nx + ny*ny
    bend[c]  = .5*np.sum(np.square(nn*plot.get_vorticity_field(nx, ny)))
    splay[c] = .5*np.sum(np.square(plot.get_gradient_field(nx, ny)))
    print "bend={}, splay={}".format(bend[c], splay[c])
    c     += 1

plt.plot(np.arange(0, c), bend)
plt.plot(np.arange(0, c), splay)
plt.show()
