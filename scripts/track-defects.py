#
# Defect tracking script. Uses python2.
#
# Usage:
#
#   python2 track_defects.py shape/nematic input output
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
from itertools import product

##################################################
# Init

if len(sys.argv) < 4:
    print("Please provide q-tensor type, input, and output files.")
    exit(1)

qtype = sys.argv[1]
if qtype not in ["shape", "nematic"]:
    print("Q-tensor type should be either 'shape' or 'nematic'.")
    exit(1)

##################################################
# Track defects


class defect_tracker:
    """
    Defect tracker.

    Honestly this is not good code... if you are reading this I feel for you!
    """

    def __init__(self, max_dst=5):
        # all defects in id order
        self.defects = []
        # list of id of active defects by charge
        self.active = {}
        # cutoff distance for matching
        self.max_dst = max_dst

    # wrapping distance on pbc domain (ugly)
    def norm(self, v1, v2, LX, LY):
        return min(
            [
                np.linalg.norm(v1 - v2 + [x * LX, y * LY])
                for x, y in product([-1, 0, 1], [-1, 0, 1])
            ]
        )

    def add_frame(self, Q00, Q01, t):
        """Add a new frame and detect defects"""
        # get dimension of pbc domain
        (LX, LY) = Q00.shape
        # get defects from new frame
        w = plot.charge_array(Q00, Q01)
        defects = plot.get_defects(w, Q00, Q01)
        # new list of active defects
        active = {}
        # compare with old ones
        for new in defects:
            # distance from all active defects with same charge
            dist = [
                self.norm(new["pos"], self.defects[i]["pos"][-1], LX, LY)
                for i in self.active.get(new["charge"], [])
            ]
            if len(dist) != 0:
                # find minimum and index
                v, i = min((v, i) for (i, v) in enumerate(dist))
                # if min dist is smaller than cutoff then we have found our defect
                if v < self.max_dst:
                    # add position and angle to stored defect
                    self.defects[self.active[new["charge"]][i]]["pos"].append(
                        new["pos"]
                    )
                    # register it as active for next round
                    active[new["charge"]] = active.get(new["charge"], []) + [
                        self.active[new["charge"]][i]
                    ]
                    # delete from current list
                    del self.active[new["charge"]][i]
                    continue

            # we have found a new defect
            self.defects.append(
                {"charge": new["charge"], "birth": t, "pos": [new["pos"]]}
            )
            # register it as active for next round
            active[new["charge"]] = active.get(new["charge"], []) + [
                len(self.defects) - 1
            ]

        # save list of acive defects for next run (all the remaining defects are inactive)
        self.active = active

    def plot(self, engine):
        # get all active defects with a certain charge
        for chg, dft in self.active.iteritems():
            # get all indices of these defects
            for ind in dft:
                p = self.defects[ind]["pos"][-1]
                engine.plot(
                    p[0],
                    p[1],
                    "b" if chg == 0.5 else "g",
                    marker=r"$ {} $".format(ind),
                    markersize=15,
                )


ar = archive.loadarchive(sys.argv[2])
tracker = defect_tracker(max_dst=15)

for i in range(0, ar._nframes + 1):
    frame = ar.read_frame(i)
    print("{}/{}".format(i, ar._nframes))
    if qtype == "shape":
        (Q00, Q01) = plot.get_nematic_field(frame.phi, frame.S00, frame.S01, size=24)
    else:
        (Q00, Q01) = plot.get_nematic_field(frame.phi, frame.Q00, frame.Q01, size=24)
    tracker.add_frame(Q00, Q01, i)

# convert format to pure np array (ok this is not great)
result = np.array(
    [[i["charge"], i["birth"], np.array(i["pos"])] for i in tracker.defects],
    dtype=object,
)
np.save(sys.argv[3], result)
