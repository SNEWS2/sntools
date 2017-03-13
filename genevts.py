#!/usr/bin/python

# call this as
# $ ./genevts.py --hierarchy [noosc|normal|inverted] --channel [ibd|es|all] -i infile -o outfile -d [SuperK|HyperK]
# where the input files are called infile_{e,eb,x}.txt and the output file is outfile.txt

from os import system
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

optchoices = ["SuperK", "HyperK"]
optdefault = "SuperK"
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
	print "channel   =", channel
	print "hierarchy =", hierarchy
	print "inputs    =", in_e, in_eb, in_x
	print "output    =", output
	print "detector  =", detector, "\n"

# call script for each interaction channel as
#     ./channel.py -i infile -o outfile -n normalization_factor -d detector
# where the normalization factor is 0 or 1 or sin^2(theta_12) or cos^2(theta_12),
# depending on oscillation scenario (see p. 236 HK public DR). We assume P_H = 0.

# normalization factors, from C. Patrignani et al. (Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
sin2t12 = 0.304
cos2t12 = 1 - sin2t12
tmpfiles = []

def execute(thisChannel, flavor, n):
	tmpfile = "tmp_%s_%s.txt" % (thisChannel, flavor)
	cmd = "python %s.py -i %s_%s.txt -o %s -n %s -d %s" % (thisChannel, input, flavor, tmpfile, n, detector)
	if verbose: print cmd
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
for filename in tmpfiles:
	f = open(filename)
	
	# read in all events
	for line in f:
		event = map(float, line.split(","))
		# `event` now is [t, pid, energy, posx, posy, posz, dirx, diry, dirz]
		events.append(event)
	
	f.close()

# sort events by first element of the list (i.e. by time)
events.sort()

outfile = open(output)
# ... and write the events to vector file (`outfile`) in this nuance-like format
outfile.close()
