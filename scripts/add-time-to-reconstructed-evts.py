"""
Add time to files with reconstructed events by comparing with original event files.
"""
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Name of input file containing 1 event per line, each consisting of `time, PID, energy, dirx, diry, dirz`. \
                                   Requires a .csv and a .reco file with same name to exist.")
args = parser.parse_args()
input = args.input
# Check that both input files exist
if not os.path.isfile(input + '.csv'):
    print "File %s.csv missing. Aborting ..." % input
    raise IOError
if not os.path.isfile(input + '.reco'):
    print "File %s.reco missing. Aborting ..." % input
    raise IOError

# Don't overwrite output file if it already exists!
if os.path.isfile(input + '.reco2'):
    print "File %s.reco2 already exists. Aborting ..." % input
    raise IOError


# Read in original (truth) events
origevents = []
with open(input + '.csv') as origfile:
    for line in origfile:
        if not (line.startswith("#") or line.isspace()):
            origevents.append(map(float, line.split(",")))

# sort them by time
origevents.sort()


# Read in reconstructed events
recoevents = []
with open(input + '.reco') as recofile:
    for line in recofile:
        if line.startswith("t,"): # ignore unnecessary hk-BONSAI output
            recoevents.append(line.split(","))

# Set time of reconstructed event to the true time. (Reconstruction uncertainty
# is at nanosecond level, negligible compared to timescale of SN simulation!)
for (i, orig) in enumerate(origevents):
    recoevents[i][0] = orig[0]
    recoevents[i][1] = 0 # PID: can't distinguish electron and positron in reconstruction


with open(input + '.reco2', 'w') as outfile:
    for (i, event) in enumerate(recoevents):
        (t, pid, ene, dirx, diry, dirz) = map(float, event)
        outfile.write("%.5f, %i, %.5f, %.5f, %.5f, %.5f\n" % (t, pid, ene, dirx, diry, dirz))
