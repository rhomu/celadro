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
    oname = "movie_"+sys.argv[2]+"_corr"
    print "Output name is", sys.argv[2]

##################################################
# plot simple animation of phases

def myplot(frame, engine):
    vx, vy = plot.get_velocity_field(frame.phi, frame.velocity, size=24)
    w      = plot.get_vorticity_field(vx, vy)
    corr1  = plot.get_corr2(vx, vy)
    corr2  = plot.get_corr(w)
    engine.plot(range(len(corr1)), corr1, label='Velocity correlation')
    engine.plot(range(len(corr2)), corr2, label='Vorticity correlation')
    engine.legend()

if len(oname)==0:
    animation.animate(ar, myplot, show=True);
else:
    an = animation.animate(ar, myplot, show=False)
    animation.save(an, oname+'.mp4', 5)
