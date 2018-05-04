#!/usr/bin/python

import __builtin__
import argparse
from datetime import datetime
from os import remove
import random

from channel import gen_evts


"""Define and parse command line options."""
parser = argparse.ArgumentParser()

parser.add_argument("input_file", help="Name or common prefix of the input file(s). Required.")

choices = ["garching", "nakazato", "totani"]
default = "totani"
parser.add_argument("-f", "--format", metavar="FORMAT", choices=choices, default=default,
                    help="Format of input files. See parsers in folder 'formats/' \
                          for details. Choices: %s. Default: %s" % (choices, default))

default = "outfile.kin"
parser.add_argument("-o", "--output", metavar="FILE", default=default,
                    help="Name of the output file. Default: '%s'." % default)

choices = ["noosc", "normal", "inverted"]
default = choices[0]
parser.add_argument("-H", "--hierarchy", metavar="HIERARCHY", choices=choices, default=default,
                    help="Oscillation scenario. Choices: %s. Default: %s" % (choices, default))

choices = ["all", "ibd", "es", "o16e", "o16eb"]
default = choices[0]
parser.add_argument("-c", "--channel", metavar="INTCHANNEL", choices=choices, default=default,
                    help="Interaction channels to consider. Currently, inverse beta decay (ibd), \
                          electron scattering (es), nu_e + oxygen CC (o16e) and nu_e-bar + oxygen CC \
                          (o16eb) are supported. Choices: %s. Default: %s" % (choices, default))

# [radius, height] of inner detector in cm
detectors = {"SuperK":[3368.15/2., 3620.],
             "HyperK":[7080./2., 5480.]}
choices = list(detectors)
default = choices[1]
parser.add_argument("-d", "--detector", metavar="DETECTOR", choices=choices, default=default,
                    help="Detector configuration. Choices: %s. Default: %s" % (choices, default))

default = 10.0
parser.add_argument("--distance", type=float, default=default,
                  help="Distance to supernova in kpc. Default: '%s'." % default)

parser.add_argument("--starttime", metavar="T", type=float,
                  help="Start generating events at T milliseconds. Default: First time bin in input file.")

parser.add_argument("--endtime", metavar="T", type=float,
                  help="Stop generating events at T milliseconds. Default: Last time bin in input file.")

parser.add_argument("-v", "--verbose", action="count",
                  help="Verbose output, e.g. for debugging. Off by default.")

args = parser.parse_args()

hierarchy = args.hierarchy
channel = args.channel
input = args.input_file
format = args.format
output = args.output
detector = args.detector
distance = args.distance
starttime = args.starttime if args.starttime else None
endtime = args.endtime if args.endtime else None
verbose = args.verbose

if verbose:
    print "channel      =", channel
    print "hierarchy    =", hierarchy
    print "input file   =", input, "--- format =", format
    print "output       =", output
    print "detector     =", detector
    print "distance     =", distance
    print "starttime    =", starttime
    print "endtime      =", endtime


'''
Event generation section.
* Handle neutrino oscillations.
* Execute scripts for all requested interaction channels.
* Interaction channel scripts generate events and write them to tmp files.

The input files (see `--input` command line option) give the fluxes of different flavors
at the boundary of the computer simulation (i.e. inside the supernova). While exiting the
supernova, neutrinos experience flavor transitions via the MSW effect.

To account for the dependence on mass hierarchy (selected via the `--hierarchy` command
line option), we run the script for each selected interaction channel on each input flux
with an appropriate normalization factor (see p. 236 of the 2016 Hyper-K Design Report).
We assume adiabatic transition (P_H = 0).
'''
# normalization factors, from C. Patrignani et al. (Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
sin2t12 = 0.304
cos2t12 = 1 - sin2t12

# List of files that contain events for a single channel/flavor combination. Will be combined into `outfile` later.
tmpfiles = []

def execute(this_channel, original_flavor, n, detected_flavor=""):

    # TODO: Replace this with a more sensible design, e.g. see https://stackoverflow.com/a/15959638
    if this_channel == "es": __builtin__._flavor = detected_flavor

    n = n * (10.0/distance)**2 # flux is proportional to 1/distance**2
    tmpfile = "tmp_%s_%s%s.txt" % (this_channel, original_flavor, detected_flavor)
    tmpfiles.append(tmpfile)

    cmd = "gen_evts(channel='%s', input='%s', format='%s', inflv='%s', output='%s', normalization=%s, detector='%s', starttime=%s, endtime=%s, verbose=%s)" \
        % (this_channel, input, format, original_flavor, tmpfile, n, detector, starttime, endtime, verbose)
    if verbose: print "Now executing:", cmd
    __builtin__._cmd = cmd
    exec(cmd)

if hierarchy == "noosc":
    if (channel == "ibd" or channel == "all"):
        execute("ibd", "eb", 1)
    if (channel == "es" or channel == "all"):
        execute("es", "e",  1, "e")
        execute("es", "eb", 1, "eb")
        execute("es", "x",  2, "x")  # normalization=2 to include both nu_mu and nu_tau
        execute("es", "xb", 2, "xb")
    if (channel == "o16e" or channel == "all"):
        execute("o16e", "e", 1)
    if (channel == "o16eb" or channel == "all"):
        execute("o16eb", "eb", 1)

if hierarchy == "normal":
    if (channel == "ibd" or channel == "all"):
        execute("ibd", "eb", cos2t12)
        execute("ibd", "xb", sin2t12)
    if (channel == "es" or channel == "all"):
        execute("es", "x",  1, "e") # nu_e that originated as nu_x
        execute("es", "eb", cos2t12, "eb") # anti-nu_e that originated as anti-nu_e
        execute("es", "xb", sin2t12, "eb") # anti-nu_e that originated as anti-nu_x
        execute("es", "e",  1, "x") # nu_x that originated as nu_e
        execute("es", "x",  1, "x") # nu_x that originated as nu_x
        execute("es", "eb", sin2t12, "xb") # anti-nu_x that originated as anti-nu_e
        execute("es", "xb", 1+cos2t12, "xb") # anti-nu_x that originated as anti-nu_x
    if (channel == "o16e" or channel == "all"):
        execute("o16e", "x", 1)
    if (channel == "o16eb" or channel == "all"):
        execute("o16eb", "eb", cos2t12)
        execute("o16eb", "xb", sin2t12)

if hierarchy == "inverted":
    if (channel == "ibd" or channel == "all"):
        execute("ibd", "xb", 1)
    if (channel == "es" or channel == "all"):
        execute("es", "e",  sin2t12, "e") # nu_e that originated as nu_e
        execute("es", "x",  cos2t12, "e") # nu_e that originated as nu_x
        execute("es", "xb", 1, "eb") # anti-nu_e that originated as anti-nu_x
        execute("es", "e",  cos2t12, "x") # nu_x that originated as nu_e
        execute("es", "x",  1+sin2t12, "x") # nu_x that originated as nu_x
        execute("es", "eb", 1, "xb") # anti-nu_x that originated as anti-nu_e
        execute("es", "xb", 1, "xb") # anti-nu_x that originated as anti-nu_x
    if (channel == "o16e" or channel == "all"):
        execute("o16e", "e", sin2t12)
        execute("o16e", "x", cos2t12)
    if (channel == "o16eb" or channel == "all"):
        execute("o16eb", "xb", 1)


'''
Output section.
* Read events generated by all interaction channel scripts.
* Distribute them randomly in the detector volume.
* Write NUANCE formatted events into output file.
'''
events = []
# read in all events:
for filename in tmpfiles:
    with open(filename) as f:
        for line in f:
            if line.startswith("#"): continue # ignore comments
            event = map(float, line.split(",")) # list(map(float, line.split(","))) in python3
            events.append(event)

    # only keep tmpfiles around in verbose mode, e.g. for verification plots
    if not verbose:
        remove(filename)

# sort events by time (i.e. the first element of the list)
events.sort()

# Write the events out in this nuance-like format
with open(output, 'w') as outfile:
    if verbose: # write parameters to file as a comment
        outfile.write("# Generated on %s with the options:\n" % datetime.now())
        outfile.write("# " + str(args) + "\n")

    for (i, event) in enumerate(events):
        (t, pid, ene, dirx, diry, dirz) = event

        # create random vertex position inside the detector volume
        rad    = detectors[detector][0] - 20.
        height = detectors[detector][1] - 20.
        while True:
            x = random.uniform(-rad, rad)
            y = random.uniform(-rad, rad)
            if x**2 + y**2 < rad**2: break
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
