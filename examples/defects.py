#
# Defect tracking script. Uses python2.
#
# Usage:
#
#   python2 defects.py input output
#
#  where
#
#    intput -- the input file or directory
#    output -- the output file

import numpy as np
import sys
sys.path.insert(0, "../plot/")
import plot
import archive
import animation

##################################################
# Init

if len(sys.argv)<3:
    print "Please provide both input and output files."
    exit(1)

##################################################
# Track defects

ar = archive.loadarchive(sys.argv[1])
tracker = plot.defect_tracker(max_dst=15)

for i in range(0, ar._nframes+1):
    frame = ar.read_frame(i)
    print "{}/{}".format(i, ar._nframes)
    (Q00, Q01) = plot.get_Qtensor(frame.phi, frame.Q00, frame.Q01, size=24)
    tracker.add_frame(Q00, Q01, i)

# convert format to pure np array (ok this is not great)
result = np.array([ [ i['charge'], i['birth'], np.array(i['pos']) ] for i in tracker.defects ],
                  dtype=object)
np.save(sys.argv[2], result)
