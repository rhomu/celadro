#
# Computes single cells trajectories and msd. Uses python2.
#
# Usage:
#
#   python2 single-cells.py input [output]
#
#  where
#
#    intput -- the input file or directory
#    output -- (optional) if present saves plots instead of showing them

import matplotlib.pyplot as plt
import sys
import numpy as np
from itertools import product

# import local libs
sys.path.insert(0, "../plot/")
import archive

##################################################
# Init

if len(sys.argv) == 1:
    print("Please provide an input file.")
    exit(1)

# load archive from file
ar = archive.loadarchive(sys.argv[1])

oname = ""
if len(sys.argv) == 3:
    oname = sys.argv[2]
    print("Output name is", sys.argv[2])

##################################################
# compute stuff and other things
size = ar.nphases
traj = [[] for i in range(size)]


# for i in np.arange(int(0.3*ar._nframes), ar._nframes+1, step=1):
for i in np.arange(0, ar._nframes+1, step=1):

    frame = ar.read_frame(i)
    print("{}/{}".format(i, ar._nframes))

    for j in range(size):
        traj[j].append(frame.com[j])


def norm(v1, v2, LX, LY):
    return min(
            [np.linalg.norm(v1-v2+[x*LX, y*LY])
             for x, y in product([-1, 0, 1], [-1, 0, 1])]
            )


# compute distance
dist = [[norm(t[i], t[0], ar.Size[0], ar.Size[1]) for i in range(len(t))] for t in traj]

for d in dist:
    plt.plot(range(len(d)), d)
plt.show()
