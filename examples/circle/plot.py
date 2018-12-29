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
import sys

sys.path.insert(0, "../plot/")

import plot
import archive
import animation

##################################################
# Init

if len(sys.argv) == 1:
    print("Please provide an input file.")
    exit(1)

# load archive from file
ar = archive.loadarchive(sys.argv[1])

oname = ""
if len(sys.argv) == 3:
    oname = "movie_"+sys.argv[2]
    print("Output name is", sys.argv[2])


##################################################
# plot simple animation of phases


def myplot(frame, fig):

    ax = fig.add_subplot(111)

    plot.cells(frame, ax)
    colors = ['b' if t == 'A' else 'r' for t in frame.parameters['types']]
    plot.solidarea(frame, ax, colors)
    plot.shape(frame, ax)
    plot.velocity(frame, ax)
    plot.pressure_force(frame, ax)

    ax.axes.set_aspect('equal', adjustable='box')
    ax.set_xlim([0, frame.parameters['Size'][0]-1])
    ax.set_ylim([0, frame.parameters['Size'][1]-1])
    ax.axis('off')


if len(oname) == 0:
    animation.animate(ar, myplot, show=True)
    exit(0)
else:
    an = animation.animate(ar, myplot, show=False)
    animation.save(an, oname+'.mp4', 5)
