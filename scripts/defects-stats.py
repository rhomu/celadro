#
# Get defect stats. Uses python2.
#
# Usage:
#
#   python2 defects-stats.py input defects
#
#  where
#
#    intput  -- the input file or directory
#    defects -- the defects file from defects.py

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
from itertools import product
from math import log
import sys

# import local libs
sys.path.insert(0, "../plot/")
import archive

##################################################
# Init

if len(sys.argv)<3:
    print "Please provide both input and defects files."
    exit(1)

# load archive from file
ar = archive.loadarchive(sys.argv[1])

# wrapping dist
def dist(v1, v2, LX, LY):
    return min(
            [ np.linalg.norm(v1-v2+[x*LX, y*LY])
              for x, y in product([-1, 0, 1], [-1, 0, 1]) ]
            )

##################################################
# Get stats

size = ar.Size[0]*ar.Size[1]

defects = np.load(sys.argv[2])

time = np.arange(ar.nstart, ar.nsteps+ar.ninfo, ar.ninfo)
dens = np.zeros(ar._nframes+1)
char = np.zeros(ar._nframes+1)
# mean speed of plus and minus defects
velp = np.zeros(len(defects))
velm = np.zeros(len(defects))

for i in range(len(defects)):

    d = defects[i]

    #print len(d[2])

    # compute density of defects and charge density
    for j in range(d[1], d[1]+len(d[2])):
        dens[j] += 1./size
        char[j] += d[0]

    # speed (only take defects that survive long enough)
    if len(d[2])>5:

        vel = np.mean([ dist(d[2][k-1], d[2][k], ar.Size[0], ar.Size[1]) for k in range(1, len(d[2])) ])

        if d[0]>0:
            velp[i] = vel
        else:
            velm[i] = vel

# mean defect density
print "Mean defect density:", np.mean(dens)
print "Mean velcotiy of plus defects:", np.mean(velp)
print "Mean velcotiy of minus defects:", np.mean(velm)

plt.figure()
plt.subplot(211)
plt.plot(time, char/size, 'r', label='charge density')
plt.plot(time, dens, label='defect density')
plt.plot(time, ndimage.uniform_filter1d(dens, size=50), label='smoothed defect density')
plt.legend()
plt.subplot(212)
plt.hist([ log(len(d[2])) for d in defects ], label='log(defect lifetime)')
plt.legend()
plt.show()
