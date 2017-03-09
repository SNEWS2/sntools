#!/usr/bin/python

# call this as
# $ ./genevts.py --channel [ibd|es|all] --hierarchy [noosc|normal|inverted] -i infile -o outfile
# where the input files are called infile_{e,eb,x}.txt and the output file is outfile.txt

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

parser.add_option("-v", "--verbose", dest="verbose",
                  help="Verbose output, e.g. for debugging. Off by default.",
                  default=False, action="store_true")

# TODO: add option for switching between SK and HK

(options, args) = parser.parse_args()

channel = options.channel
hierarchy = options.hierarchy
in_e = options.input + "_e.txt"
in_eb = options.input + "_eb.txt"
in_x = options.input + "_x.txt"
output = options.output
verbose = options.verbose

if verbose:
	print "channel   =", channel
	print "hierarchy =", hierarchy
	print "inputs    =", in_e, in_eb, in_x
	print "output    =", out, "\n"

# call script for each interaction channel as
#     ./channel.py -i infile -o outfile -n normalization_factor
# where the normalization factor is 0 or 1 or sin^2(theta_12) or cos^2(theta_12),
# depending on oscillation scenario (see p. 236 HK public DR). We assume P_H = 0.

# normalization factors, from C. Patrignani et al. (Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
sin2t12 = 0.304
cos2t12 = 1 - sin2t12
tmpfiles = []

if (hierarchy == "noosc"):
	if (channel == "ibd" or channel == "all"):
		print "./ibd.py -i %s -o tmp_eb_ibd -n 1" % in_eb
		tmpfiles.append("tmp_eb_ibd")
	if (channel == "es" or channel == "all"):
		print "./es.py -i %s -o tmp_e_es -n 1" % in_e
		tmpfiles.append("tmp_e_es")

if (hierarchy == "normal"):
	if (channel == "ibd" or channel == "all"):
		print "./ibd.py -i %s -o tmp_eb_ibd -n cos2t12" % in_eb
		print "./ibd.py -i %s -o tmp_x_ibd -n sin2t12" % in_x
		tmpfiles.extend(["tmp_eb_ibd","tmp_x_ibd"])
	if (channel == "es" or channel == "all"):
		print "./es.py -i %s -o tmp_x_es -n 1" % in_x
		tmpfiles.append("tmp_x_es")

if (hierarchy == "inverted"):
	if (channel == "ibd" or channel == "all"):
		print "./ibd.py -i %s -o tmp_x_ibd -n 1" % in_x
		tmpfiles.append("tmp_x_ibd")
	if (channel == "es" or channel == "all"):
		print "./es.py -i %s -o tmp_e_es -n sin2t12" % in_e
		print "./es.py -i %s -o tmp_x_es -n cos2t12" % in_x
		tmpfiles.extend(["tmp_e_es","tmp_x_es"])

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
