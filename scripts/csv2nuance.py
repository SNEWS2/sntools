"""
Batch convert csv event files to NUANCE-formatted files that can be
used as input files to WCSim.
Can generate a WCSim macro file that simulates these NUANCE files.
"""
import argparse
import os
import random

from genevts import detectors

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Name of input file containing 1 event per line, each consisting of `time, PID, energy, dirx, diry, dirz`. \
                                   If `input` is a directory, iterate over all `.csv` and `.txt` files in the directory.")
parser.add_argument("-o", "--output", help="Name of output file in NUANCE format. Ignored if `input` is a directory. \
                                            Default: same as input file, but with .kin extension.")
parser.add_argument("-w", "--wcsim", help="Write WCSim macro file containing all converted files.")
choices = list(detectors)
default = choices[1]
parser.add_argument("-d", "--detector", metavar="DETECTOR", choices=choices, default=default,
                    help="Detector configuration. Choices: %s. Default: %s" % (choices, default))
args = parser.parse_args()

if os.path.isfile(args.input):
    inputs = [args.input]
elif os.path.isdir(args.input):
    inputs = [os.path.join(args.input, file) for file in os.listdir(args.input)
                                             if os.path.splitext(file)[1] in (".csv", ".txt")]
else:
    raise IOError

if args.wcsim:
    wcsimfile = open(args.wcsim, 'w')
    wcsimfile.write('# configuration options (from the sample macro file in WCSim 1.7.0)\n')
    wcsimfile.write('/run/verbose 0\n')
    wcsimfile.write('/tracking/verbose 0\n')
    wcsimfile.write('/hits/verbose 0\n')
    wcsimfile.write('/WCSim/WCgeom ' + args.detector + '\n')
    wcsimfile.write('/WCSim/Construct\n')
    wcsimfile.write('/WCSim/PMTQEMethod Stacking_Only\n')
    wcsimfile.write('/WCSim/PMTCollEff on\n')
    wcsimfile.write('/WCSim/SavePi0 false\n')
    wcsimfile.write('/DAQ/Digitizer SKI\n')
    wcsimfile.write('/DAQ/Trigger NDigits\n')
    wcsimfile.write('/control/execute daq.mac\n')
    wcsimfile.write('/DarkRate/SetDarkMode 1\n')
    wcsimfile.write('/DarkRate/SetDarkHigh 100000\n')
    wcsimfile.write('/DarkRate/SetDarkLow 0\n')
    wcsimfile.write('/DarkRate/SetDarkWindow 4000\n')
    wcsimfile.write('/mygen/generator muline\n\n')
    wcsimfile.write('# Vectorfiles to simulate:\n')

for input in inputs:
    # determine name of output file
    if args.output and os.path.isfile(args.input):
        output = args.output
    else:
        root, ext = os.path.splitext(input)
        output = root + ".kin"

    # make sure we don't overwrite an existing file
    if os.path.isfile(output):
        print "File %s already exists. Aborting ..." % output
        exit()

    # read from input file
    events, comments = [], []
    with open(input) as infile:
        for line in infile:
            if (line.startswith("#") or line.isspace()):
                comments.append(line)
            else:
                events.append(map(float, line.split(",")))

    events.sort() # sort by time (i.e. the first element of the list)

    # write events to output file
    with open(output, 'w') as outfile:
        for line in comments:
            # Need to wait until WCSim ignores comment lines.
            # See https://github.com/WCSim/WCSim/issues/234
            pass #outfile.write(line)

        for (i, event) in enumerate(events):
            (t, pid, ene, dirx, diry, dirz) = event

            # create random vertex position inside the detector volume
            radius = detectors[args.detector][0] - 20
            height = detectors[args.detector][1] - 20
            while True:
                x = random.uniform(-radius, radius)
                y = random.uniform(-radius, radius)
                if x**2 + y**2 < radius**2: break
            z = random.uniform(-height/2, height/2)

            outfile.write("$ begin\n")
            outfile.write("$ nuance 0\n")
            outfile.write("$ vertex %.5f %.5f %.5f %.5f\n" % (x, y, z, t))
            outfile.write("$ track 14 1020.00000 1.00000 0.00000 0.00000 -1\n") # "Neutrino" Track
            outfile.write("$ track 2212 935.98400 0.00000 0.00000 1.00000 -1\n") # "Target" track
            outfile.write("$ info 0 0 %i\n" % i)
            outfile.write("$ track %i %.5f %.5f %.5f %.5f 0\n" % (pid, ene, dirx, diry, dirz)) # Outgoing particle track
            outfile.write("$ end\n")

        outfile.write("$ stop\n")

    # write WCSim macro file
    if args.wcsim:
        wcsimfile.write('/mygen/vecfile ' + output + '\n')
        wcsimfile.write('/WCSimIO/RootFile ' + os.path.splitext(input)[0] + '.root' + '\n')
        wcsimfile.write('/run/beamOn ' + str(len(events)) + '\n\n')

if args.wcsim:
    wcsimfile.close()
