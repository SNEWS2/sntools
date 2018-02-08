#!/usr/bin/python

# call this as
# $ ./genevts.py --hierarchy [noosc|normal|inverted] --channel [ibd|es|all] -i infile -o outfile.kin -d [SuperK|HyperK]
# where the input files are called infile_{e,eb,x}.txt and the output file is outfile.kin

import random
from optparse import OptionParser
import __builtin__
import channel as chnl


'''
Setup section.
* Define command line options.
* Parse input options (and perform some sanity checks).
'''
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
optchoices = detectors.keys() # list(detectors.keys()) in python3
optdefault = detectors.keys()[0]
parser.add_option("-d", "--detector", dest="detector",
                  help="Detector configuration. Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
                  metavar="DETECTOR",
                  choices=optchoices, default=optdefault)

parser.add_option("--starttime", dest="starttime",
                  help="Start generating events at T seconds. Useful to speed up calculation if you are only interested in a short time window. Default: First time bin in input file.",
                  metavar="T")

parser.add_option("--endtime", dest="endtime",
                  help="Stop generating events at T seconds. Useful to speed up calculation if you are only interested in a short time window. Default: Last time bin in input file.",
                  metavar="T")

parser.add_option("-v", "--verbose", dest="verbose",
                  help="Verbose output, e.g. for debugging. Off by default.",
                  default=False, action="store_true")

(options, args) = parser.parse_args()

hierarchy = options.hierarchy
channel = options.channel
input = options.input
output = options.output
detector = options.detector
starttime = options.starttime if options.starttime else None
endtime = options.endtime if options.endtime else None
verbose = options.verbose

if verbose:
	print "channel   =", channel
	print "hierarchy =", hierarchy
	print "inputs    =", input + "_e.txt", input + "_eb.txt", input + "_x.txt"
	print "output    =", output
	print "detector  =", detector
	print "starttime =", starttime
	print "endtime   =", endtime


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
	
	infile = "%s_%s.txt" % (input, original_flavor)
	tmpfile = "tmp_%s_%s%s.txt" % (this_channel, original_flavor, detected_flavor)
	tmpfiles.append(tmpfile)
	
	cmd = "chnl.main(channel='%s', input='%s', output='%s', normalization=%s, detector='%s', starttime=%s, endtime=%s, verbose=%s)" % (this_channel, infile, tmpfile, n, detector, starttime, endtime, verbose)
	if verbose: print "Now executing:", cmd
	exec(cmd)

if (hierarchy == "noosc"):
	if (channel == "ibd" or channel == "all"):
		execute("ibd", "eb", 1)
	if (channel == "es" or channel == "all"):
		execute("es", "e",  1, "e")
		execute("es", "eb", 1, "eb")
		execute("es", "x",  2, "x")  # normalization=2 to include both nu_mu and nu_tau
		execute("es", "x",  2, "xb") # anti-nu_x have different cross section then nu_x but use same input file

if (hierarchy == "normal"):
	if (channel == "ibd" or channel == "all"):
		execute("ibd", "eb", cos2t12)
		execute("ibd", "x",  sin2t12)
	if (channel == "es" or channel == "all"):
		execute("es", "x",  1, "e") # nu_e that originated as nu_x
		execute("es", "eb", cos2t12, "eb") # anti-nu_e that originated as anti-nu_e
		execute("es", "x",  sin2t12, "eb") # anti-nu_e that originated as anti-nu_x
		execute("es", "e",  1, "x") # nu_x that originated as nu_e
		execute("es", "x",  1, "x") # nu_x that originated as nu_x
		execute("es", "eb", sin2t12, "xb") # anti-nu_x that originated as anti-nu_e
		execute("es", "x",  1+cos2t12, "xb") # anti-nu_x that originated as anti-nu_x

if (hierarchy == "inverted"):
	if (channel == "ibd" or channel == "all"):
		execute("ibd", "x", 1)
	if (channel == "es" or channel == "all"):
		execute("es", "e",  sin2t12, "e") # nu_e that originated as nu_e
		execute("es", "x",  cos2t12, "e") # nu_e that originated as nu_x
		execute("es", "x",  1, "eb") # anti-nu_e that originated as anti-nu_x
		execute("es", "e",  cos2t12, "x") # nu_x that originated as nu_e
		execute("es", "x",  1+sin2t12, "x") # nu_x that originated as nu_x
		execute("es", "eb", 1, "xb") # anti-nu_x that originated as anti-nu_e
		execute("es", "x",  1, "xb") # anti-nu_x that originated as anti-nu_x


'''
Output section.
* Read events generated by all interaction channel scripts.
* Distribute them randomly in the detector volume.
* Write NUANCE formatted events into output file.
'''
events = [] # this will become a list of lists: one entry per event, which is a list of time, energy, etc.
# read in all events:
for filename in tmpfiles:
	f = open(filename)
	for line in f:
		event = map(float, line.split(",")) # list(map(float, line.split(","))) in python3
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
	
	if verbose: print "events[%d] = %s" %(i, event)
	
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
