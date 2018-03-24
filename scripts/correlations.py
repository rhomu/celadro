#
# Computes mean correlations (vorticity^2, velocity^2, Q^2). Uses python2.
#
# Usage:
#
#   python2 correlations.py input [output]
#
#  where
#
#    intput -- the input file or directory
#    output -- (optional) if present saves plots instead of showing them

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

oname = ""
if len(sys.argv)==3:
    oname = sys.argv[2]
    print "Output name is", sys.argv[2]

##################################################
# compute mean correlations

size = int(sqrt(ar.Size[0]*ar.Size[1])/2)
count = 0
vel_tot = np.zeros(size)
vel_err = np.zeros(size)
vor_tot = np.zeros(size)
vor_err = np.zeros(size)
nem_tot = np.zeros(size)
nem_err = np.zeros(size)

for i in np.arange(int(0.3*ar._nframes), ar._nframes+1, step=1):

    frame = ar.read_frame(i)
    print "{}/{}".format(i, ar._nframes)

    # velocity and vorticity
    vx, vy = plot.get_velocity_field(frame.phi, frame.velocity, size=24)
    w      = plot.get_vorticity_field(vx, vy)
    corr1  = plot.get_corr2(vx, vy)
    corr2  = plot.get_corr(w)

    vel_tot += corr1
    vel_err += np.square(corr1)
    vor_tot += corr2
    vor_err += np.square(corr2)

    # Q
    Qxx, Qxy = plot.get_Qtensor(frame.phi, frame.Q00, frame.Q01, size=24)
    corrQ    = plot.get_corr2(Qxx, Qxy) # nornalization is wrong but cancels

    nem_tot += corrQ
    nem_err += np.square(corrQ)

    count += 1

for i in range(size):
    vel_tot[i] /= count
    vel_err[i] /= count
    vel_err[i]  = sqrt(vel_err[i] - vel_tot[i]**2)
    vor_tot[i] /= count
    vor_err[i] /= count
    vor_err[i]  = sqrt(vor_err[i] - vor_tot[i]**2)
    nem_tot[i] /= count
    nem_err[i] /= count
    nem_err[i]  = sqrt(nem_err[i] - nem_tot[i]**2)

#
# nematic
#
plt.figure()
plt.plot(range(size), nem_tot, 'k')
plt.fill_between(range(size), nem_tot-nem_err, nem_tot+nem_err)
if oname=='': plt.show()
else:
    plt.savefig(oname+'_corr_nem.png')
    np.save(oname+'_nem_tot', nem_tot)
    np.save(oname+'_nem_err', nem_err)

#
# velocity
#
plt.figure()
plt.plot(range(size), vel_tot, 'k')
plt.fill_between(range(size), vel_tot-vel_err, vel_tot+vel_err)
if oname=='': plt.show()
else:
    plt.savefig(oname+'_corr_vel.png')
    np.save(oname+'_vel_tot', vel_tot)
    np.save(oname+'_vel_err', vel_err)

#
# vorticity
#
def zero(x, y):
    indi = np.where(y[1:]*y[0:-1] < 0.0)[0][0]
    dx = x[indi+1] - x[indi]
    dy = y[indi+1] - y[indi]
    z  = -y[indi] * (dx/dy) + x[indi]
    return z

# get location of zero of the derivative
dvor_tot = np.gradient(vor_tot)
xz = zero(range(size), dvor_tot)
yz = vor_tot[int(xz)]+(xz-int(xz))*(vor_tot[int(xz)+1]-vor_tot[int(xz)])

plt.figure()
plt.plot(range(size), vor_tot, 'k')
plt.plot(xz, yz, 'ro')
plt.fill_between(range(size), vor_tot-vor_err, vor_tot+vor_err)
if oname=='': plt.show()
else:
    plt.savefig(oname+'_corr_vor.png')
    np.save(oname+'_vor_tot', vor_tot)
    np.save(oname+'_vor_err', vor_err)
