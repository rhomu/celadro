#
# This is a simple example file to show the plotting capabilities of the
# program.
#
# Usage:
#
#   python plot.py input [output]
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
    oname = "movie_" + sys.argv[2]
    print("Output name is", sys.argv[2])


##################################################
# plot simple animation of phases


def myplot(frame, fig):

    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    plot.cells(frame, ax1)
    plot.nematic(frame, ax1)

    plot.nematic_field(frame, engine=ax2, avg=2, show_def=True)

    plot.substrate_phi(frame, engine=ax3, cbar=True)

    for ax in [ax1, ax2, ax3]:
        ax.axes.set_aspect("equal", adjustable="box")
        ax.set_xlim([0, frame.parameters["Size"][0] - 1])
        ax.set_ylim([0, frame.parameters["Size"][1] - 1])
        ax.axis("off")


if len(oname) == 0:
    animation.animate(ar, myplot, show=True)
else:
    an = animation.animate(ar, myplot, show=False)
    animation.save(an, oname + ".mp4", 5)
