#!/usr/bin/python

# call this as
# $ ./genevts.py --hierarchy [noosc|normal|inverted] --channel [ibd|es|all] -i infile -o outfile -d [SuperK|HyperK]
# where the input files are called infile_{e,eb,x}.txt and the output file is outfile.txt

from os import system
import random
from optparse import OptionParser
parser = OptionParser()

optchoices = ["noosc", "normal", "inverted"]
optdefault = "noosc"
parser.add_option("-H", "--hierarchy", dest="hierarchy",
                  help="Oscillation scenario to consider. Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
                  metavar="HIERARCHY",
                  choices=optchoices, default=optdefault)

optchoices = ["all", "ibd", "es"]
optdefault = "all"
parser.add_option("-c", "--channel", dest="channel",
                  help="Interaction channels to consider. Currently, inverse beta decay (ibd) and electron scattering (es) are supported. Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
                  metavar="INTCHANNEL",
                  choices=optchoices, default=optdefault)

optdefault = "infile"
parser.add_option("-i", "--input", dest="input",
                  help="Common prefix of the input files. Default: '%s', which results in the files '%s_e.txt', '%s_eb.txt' and '%s_x.txt' being read." \
                      % (optdefault, optdefault, optdefault, optdefault),
                  metavar="PREFIX",
                  default=optdefault)

optdefault = "outfile.kin"
parser.add_option("-o", "--output", dest="output",
                  help="Name of the output file. Default: '%s'." \
                      % (optdefault),
                  metavar="FILENAME",
                  default=optdefault)

# [radius, height] of inner detector in cm
detectors = {"SuperK":[3368.15/2., 3620.],
             "HyperK":[7080./2., 5480.]}
optchoices = detectors.keys()
optdefault = detectors.keys()[0]
parser.add_option("-d", "--detector", dest="detector",
                  help="Detector configuration. Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
                  metavar="DETECTOR",
                  choices=optchoices, default=optdefault)

parser.add_option("-v", "--verbose", dest="verbose",
                  help="Verbose output, e.g. for debugging. Off by default.",
                  default=False, action="store_true")

(options, args) = parser.parse_args()

hierarchy = options.hierarchy
channel = options.channel
input = options.input
in_e = input + "_e.txt"
in_eb = input + "_eb.txt"
in_x = input + "_x.txt"
output = options.output
detector = options.detector
verbose = options.verbose

if verbose:
	print ("channel   =", channel)
	print ("hierarchy =", hierarchy)
	print ("inputs    =", in_e, in_eb, in_x)
	print ("output    =", output)
	print ("detector  =", detector, "\n")

# call script for each interaction channel as
#     ./channel.py -i infile -o outfile -n normalization_factor -d detector
# where the normalization factor is 0 or 1 or sin^2(theta_12) or cos^2(theta_12),
# depending on oscillation scenario (see p. 236 HK public DR). We assume P_H = 0.

# normalization factors, from C. Patrignani et al. (Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
sin2t12 = 0.304
cos2t12 = 1 - sin2t12

# List of files that contain events for a single channel/flavor combination. Will be combined into `outfile` later.
tmpfiles = []

def execute(thisChannel, flavor, n):
	tmpfile = "tmp_%s_%s.txt" % (thisChannel, flavor)
	cmd = "python %s.py --input=%s_%s.txt --output=%s --normalization=%s --detector=%s" % (thisChannel, input, flavor, tmpfile, n, detector)
	if verbose:
		cmd = cmd + " --verbose" # inherit verbosity
		print ("Now executing:", cmd)
	system(cmd)
	tmpfiles.append(tmpfile)

if (hierarchy == "noosc"):
	if (channel == "ibd" or channel == "all"):
		execute("ibd", "eb", 1)
	if (channel == "es" or channel == "all"):
		execute("es", "e", 1)

if (hierarchy == "normal"):
	if (channel == "ibd" or channel == "all"):
		execute("ibd", "eb", cos2t12)
		execute("ibd", "x", sin2t12)
	if (channel == "es" or channel == "all"):
		execute("es", "x", 1)

if (hierarchy == "inverted"):
	if (channel == "ibd" or channel == "all"):
		execute("ibd", "x", 1)
	if (channel == "es" or channel == "all"):
		execute("es", "e", sin2t12)
		execute("es", "x", cos2t12)

events = [] # this will become a list of lists: one entry per event, which is a list of time, energy, etc.
# read in all events:
for filename in tmpfiles:
	f = open(filename)
	for line in f:
		event = map(float, line.split(","))
		# `event` has the format `[t, pid, energy, dirx, diry, dirz]`
		events.append(event)
	f.close()

# sort events by first element of the list (i.e. by time)
events.sort()

# ... and write the events to vector file (`outfile`) in this nuance-like format
outfile = open(output, 'w')
for i in range(len(events)):
	event = events[i]
	t   = event[0]
	pid = int(event[1])
	ene = event[2]
	(dirx, diry, dirz) = (event[3], event[4], event[5])
	
	if verbose: print ("events[",i,"] = ", event)
	
	# create random vertex position inside the detector volume
	rad    = detectors[options.detector][0] - 20.
	height = detectors[options.detector][1] - 20.
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
outfile.close()
