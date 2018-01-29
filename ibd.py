#!/usr/bin/python

from optparse import OptionParser
from math import pi, sin, cos, sqrt, gamma, exp, floor, ceil
from scipy import integrate, interpolate
import numpy as np

'''
Setup section.
* Define command line options.
* Parse input options (and perform some sanity checks).
* Read in data from input file.
'''
parser = OptionParser()

optdefault = "infile_eb.txt"
parser.add_option("-i", "--input", dest="input",
                  help="Name of the input file. Default: '%s'." \
                      % (optdefault),
                  metavar="FILENAME",
                  default=optdefault)

optdefault = "tmp_ibd_eb.txt"
parser.add_option("-o", "--output", dest="output",
                  help="Name of the output file. Default: '%s'." \
                      % (optdefault),
                  metavar="FILENAME",
                  default=optdefault)

optdefault = 1.0
parser.add_option("-n", "--normalization", dest="normalization",
                  help="Normalization factor to account for neutrino oscillations. Gets set by `genevts.py`. Default: '%s'." \
                      % (optdefault),
                  metavar="NORMALIZATION",
                  default=optdefault)

# number of free protons (i.e. H nuclei) in each detector
detectors = {"SuperK": 2.1e+33,
             "HyperK": 1.4e+34} # one-tank configuration
optchoices = detectors.keys() # list(detectors.keys()) in python3
optdefault = detectors.keys()[0]
parser.add_option("-d", "--detector", dest="detector",
                  help="Detector configuration. Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
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

verbose = options.verbose
starttime = options.starttime
endtime = options.endtime
normalization = float(options.normalization)
if (normalization <= 0 or normalization > 1):
	print("Error: Normalization factor should be in the interval (0,1]. Aborting ...")
	exit()

# read data from input file, remove lines with comments and empty lines
with open(options.input) as infile:
	if verbose: print "Reading neutrino simulation data from", options.input, "..."
	indata = [line for line in infile if not (line.startswith("#") or line.isspace())]

# if start time and end time are not given as command line arguments, get them from 1st/last line of input file
_starttime = indata[0].split(",")[0]
_endtime = indata[-1].split(",")[0]

if not starttime:
	starttime = _starttime
elif starttime < _starttime:
	print("Error: Start time must be greater than time in first line of input file. Aborting ...")
	exit()

if not endtime:
	endtime = _endtime
elif endtime > _endtime:
	print("Error: End time must be less than time in last line of input file. Aborting ...")
	exit()

starttime = ceil(float(starttime) * 1000) # convert to ms
endtime = floor(float(endtime) * 1000)
duration = endtime - starttime


'''
Particle physics section.
* cross section for neutrino-electron scattering
* directionality of scattered electron

Based on Strumia/Vissani (2003), arXiv:astro-ph/0302055.
'''
mN = 939.5654 # neutron mass (MeV)
mP = 938.2721 # proton mass (MeV)
mE = 0.5109989 # electron mass (MeV)
mPi = 139.57018 # pion mass (MeV)
delta = mN-mP
mAvg=(mP+mN)/2
gF=1.16637e-11 # Fermi coupling constant

def t_eNu_eE(eNu, eE):
	return mN**2 - mP**2 - 2*mP*(eNu-eE)
def x(eNu, eE):
	return t_eNu_eE(eNu, eE) / (4*mAvg**2)
def y(eNu, eE):
	return 1 - t_eNu_eE(eNu, eE)/710000
def z(eNu, eE):
	return 1 - t_eNu_eE(eNu, eE)/1000000
def f1(eNu, eE):
	return (1 - 4.706*x(eNu, eE)) / ((1-x(eNu, eE)) * y(eNu, eE)**2)
def f2(eNu, eE):
	return 3.706 / ((1-x(eNu, eE)) * y(eNu, eE)**2)
def g1(eNu, eE):
	return -1.27 / z(eNu, eE)**2
def g2(eNu, eE):
	return 2 * g1(eNu, eE) * mAvg**2 / (mPi**2 - t_eNu_eE(eNu, eE)) 

# AM, BM and CM for approximate calculation of absMsquared,
# AM1, BM1 and CM1 for more precise calculation
def AM(eNu, eE):
	return (mAvg**2 * (f1(eNu, eE)**2 - g1(eNu, eE)**2) * (t_eNu_eE(eNu, eE)-mE**2)) - (mAvg**2 * delta**2 * (f1(eNu, eE)**2 + g1(eNu, eE)**2)) - (2 * mE**2 * mAvg * delta * g1(eNu, eE) *(f1(eNu, eE)+f2(eNu, eE)))
def AM1(eNu, eE):
	return  1./16 * ( 
	(t_eNu_eE(eNu, eE) - mE**2) * (
		4 * f1(eNu, eE)**2 * (4*mAvg**2 + t_eNu_eE(eNu, eE) + mE**2)
		+ 4 * g1(eNu, eE)**2 * (-4*mAvg**2 + t_eNu_eE(eNu, eE) + mE**2)
		+ f2(eNu, eE)**2 * ((t_eNu_eE(eNu, eE)**2)/(mAvg**2) + 4*t_eNu_eE(eNu, eE) + 4*mE**2)
		+ 4*mE**2 * t_eNu_eE(eNu, eE) * g2(eNu, eE)**2 / mAvg**2
		+ 8*f1(eNu, eE)*f2(eNu, eE) * (2*t_eNu_eE(eNu, eE) + mE**2)
		+ 16*mE**2 * g1(eNu, eE)*g2(eNu, eE))
	- delta**2 * (
		(4*f1(eNu, eE)**2 + t_eNu_eE(eNu, eE) * f2(eNu, eE)**2 / mAvg**2) *
		(4*mAvg**2 + t_eNu_eE(eNu, eE) - mE**2)
		+ 4*g1(eNu, eE)**2 * (4*mAvg**2 - t_eNu_eE(eNu, eE) + mE**2)
		+ 4*mE**2 * g2(eNu, eE)**2 * (t_eNu_eE(eNu, eE) - mE**2) / mAvg**2
		+ 8*f1(eNu, eE)*f2(eNu, eE) * (2*t_eNu_eE(eNu, eE) - mE**2)
		+ 16*mE**2 * g1(eNu, eE)*g2(eNu, eE))
	- 32*mE**2 * mAvg * delta * g1(eNu, eE)*(f1(eNu, eE) + f2(eNu, eE)))

def BM(eNu, eE):
	return t_eNu_eE(eNu, eE)*g1(eNu, eE)*(f1(eNu, eE)+f2(eNu, eE))
def BM1(eNu, eE):
	return 1./16 * (
	16*t_eNu_eE(eNu, eE) * g1(eNu, eE)*(f1(eNu, eE) + f2(eNu, eE))
	+ 4*mE**2 * delta * (f2(eNu, eE)**2 + f1(eNu, eE)*f2(eNu, eE) + 2*g1(eNu, eE)*g2(eNu, eE))/mAvg)

def CM(eNu, eE):
	return ((f1(eNu, eE)**2) + (g1(eNu, eE)**2))/4
def CM1(eNu, eE):
	return 1./16 * (4*(f1(eNu, eE)**2 + g1(eNu, eE)**2) - t_eNu_eE(eNu, eE) * f2(eNu, eE)**2 / mAvg**2)

def sMinusU(eNu, eE):
	return 2*mP*(eNu+eE) - mE**2

def absMsquared(eNu, eE):
	return AM(eNu, eE) - sMinusU(eNu, eE) * BM(eNu, eE) + sMinusU(eNu, eE)**2 * CM(eNu, eE)

def dSigmadE(eNu, eE):
	return 2 * mP * gF**2 * 0.9746**2 / (8 * pi * mP**2 * eNu**2) * absMsquared(eNu, eE)


# return direction of outgoing positron, if incoming neutrino moves in z direction
def direction(eneNu):
	pMax = 0
	cosT = 0
	nCosTBins = 1000
	cosTBinWidth = 2./nCosTBins
	for j in range(nCosTBins):
		cosT = -1 + cosTBinWidth*(j+0.5) # 1000 steps in the interval [-1,1]
		p = dir_nuebar_p_sv(eneNu, cosT)
		if p > pMax:
			pMax = p
	
	while (True):
		cosT = 2*np.random.random() - 1 # randomly distributed in interval [-1,1)
		if dir_nuebar_p_sv(eneNu, cosT) > pMax*np.random.random():
			sinT = sin(np.arccos(cosT))
			phi = 2 * pi * np.random.random() - pi # randomly distributed in [-pi, pi)
			break
	
	return (sinT*cos(phi), sinT*sin(phi), cosT)

# probability distribution for the angle at which the positron is emitted
# numerical values are from Ishino-san's code for SK, based on email from Vissani
def dir_nuebar_p_sv(eneNu, cosT):
	def dir_f1(eneNu):
		return -0.05396 + 0.35824 * (eneNu/100) + 0.03309 * (eneNu/100)**2
	def dir_f2(eneNu):
		return  0.00050 - 0.02390 * (eneNu/100) + 0.14537 * (eneNu/100)**2
	
	return 0.5 + dir_f1(eneNu) * cosT + dir_f2(eneNu) * (cosT**2 -1./3)


'''
Astrophysics section.
* neutrinos are well described by a Gamma distribution
* calculate energy-dependent flux at a fiducial distance of 10 kpc
'''
# Convert fiducial distance of 10 kpc into units of MeV**(-1)
# see http://www.wolframalpha.com/input/?i=10+kpc%2F(hbar+*+c)+in+MeV%5E(-1)
dSquared = (1.563738e+33)**2 

# energy distribution of neutrinos
def gamma_dist(eNu):
	return eNu**alpha / gamma(alpha + 1) * ((alpha + 1)/a)**(alpha + 1) * exp(-(alpha + 1)* eNu/a)

def dFluxdE(eNu, luminosity):
	return 1/(4*pi*dSquared) * luminosity/a * gamma_dist(eNu)


'''
Preparation section.
* Parse input data.
* For each time step in the input data, calculate instantaneous event rate.
* Interpolate to get event rate as a function of time.
'''
# double differential event rate
def ddEventRate(eE, eNu):
	return dSigmadE(eNu, eE)*dFluxdE(eNu, L)

# calculate range for eE from eNu in center-of-mass (cm) frame
def s(eNu):
	return 2*mP*eNu + mP**2
def pE_cm(eNu):
	return (sqrt((s(eNu)-(mN-mE)**2)*(s(eNu)-(mN+mE)**2)))/(2*sqrt(s(eNu)))
def eE_cm(eNu):
	return (s(eNu)-mN**2+mE**2)/(2*sqrt(s(eNu)))

# Bounds for integration over eE and eNu
delta_cm = (mN**2 - mP**2 - mE**2)/(2*mP)
eThr=((mN+mE)**2 - mP**2)/(2*mP) # threshold energy for IBD

def eE_Min(eNu):
	return eNu - delta_cm - (eNu/sqrt(s(eNu)) * (eE_cm(eNu) + pE_cm(eNu)))
def eE_Max(eNu):
	return eNu - delta_cm - (eNu/sqrt(s(eNu)) *(eE_cm(eNu) - pE_cm(eNu)))

def bounds_eE(eNu):
	return [eE_Min(eNu)+1, eE_Max(eNu)+1]
def bounds_eNu():
	return [eThr,100]

nevtValues=[]
tValues=[]
aValues=[]
eNuSquaredValues=[]
totnevt = 0
nP = detectors[options.detector] # number of protons in detector

for line in indata:
	# get time, mean energy, mean squared energy, luminosity
	t, a, eNuSquared, L = line.split(",")
	t=float(t) * 1000 # convert to ms
	tValues.append(t)
	a=float(a)
	aValues.append(a)
	eNuSquared = float(eNuSquared)
	eNuSquaredValues.append(eNuSquared)
	L=float(L) * 624.151 # convert from erg/s to MeV/ms
	
	alpha = (2*a**2 - eNuSquared) / (eNuSquared - a**2)
	
	# integrate over eE and then eNu to obtain the event rate at time t
	simnevt = nP * integrate.nquad(ddEventRate, [bounds_eE, bounds_eNu]) [0]
	
	# create a list of event rates at time t for input into interpolation function
	nevtValues.append(simnevt)

# interpolate the mean energy, mean squared energy and event rate
interpolatedEnergy = interpolate.pchip(tValues, aValues)
interpolatedMSEnergy = interpolate.pchip(tValues, eNuSquaredValues)
interpolatedNevt = interpolate.pchip(tValues, nevtValues)


'''
Event generation section.
* For each time bin, get number of events from a Poisson distribution.
* Generate random events with appropriate distribution of time/energy/direction.
* Write them to output file.
'''
binWidth = 1 # bin width in ms
binNr = np.arange(1, floor(duration/binWidth)+1) # number of full-width bins
if verbose:
	print "Now generating events in %s ms bins between %s-%s ms" % (binWidth, starttime, endtime)
	print "**************************************"

outfile = open(options.output, 'w')
# integrate event rate and energy over each bin
for i in binNr:
	boundsMin = starttime + (i-1)*binWidth
	boundsMax = starttime + i*binWidth
	
	# calculate expected number of events in this bin and multiply with a normalization factor
	# (1, sin^2(theta_12), cos^2(theta_12)) to take neutrino oscillations into account
	binnedNevt = integrate.quad(interpolatedNevt, boundsMin, boundsMax)[0] * normalization
	# randomly select number of events in this bin from Poisson distribution around binnedNevt:
	binnedNevtRnd = np.random.choice(np.random.poisson(binnedNevt, size=1000))
	# find the total number of events over all bins
	totnevt += binnedNevtRnd
	
	# create binned values for energy, mean squared energy and shape parameter
	binnedEnergy = integrate.quad(interpolatedEnergy, boundsMin, boundsMax)[0] / binWidth
	binnedMSEnergy = integrate.quad(interpolatedMSEnergy, boundsMin, boundsMax)[0] / binWidth
	binnedAlpha = (2*binnedEnergy**2 - binnedMSEnergy)/(binnedMSEnergy - binnedEnergy**2)
	
	if verbose:
		print "timebin       = %s-%s ms" % (boundsMin, boundsMax)
		print "Nevt (theor.) =", binnedNevt
		print "Nevt (actual) =", binnedNevtRnd
		print "mean energy   =", binnedEnergy, "MeV"
		print "**************************************"
	
	# define particle for each event in time interval
	for i in range(binnedNevtRnd):
		t = boundsMin + np.random.random() * binWidth
		
		# generate a neutrino energy above eThr
		while (True):
			eneNu = np.random.gamma(binnedAlpha + 1, binnedEnergy/(binnedAlpha + 1))
			if eneNu > eThr:
				break
		
		# generate direction of positron at given neutrino energy
		(dirx, diry, dirz) = direction(eneNu)
		# generate positron energy at given neutrino energy and cosT (Strumia & Vissani, 2003)
		epsilon = eneNu/mP
		kappa = (1 + epsilon)**2 - (epsilon * dirz)**2
		ene = ((eneNu - delta_cm) * (1 + epsilon) + epsilon * dirz * sqrt((eneNu - delta_cm)**2 - mE**2 * kappa)) / kappa
		# print out [t, PID, energy, dirx, diry, dirz] to file
		outfile.write("%f, -11, %f, %f, %f, %f\n" % (t, ene, dirx, diry, dirz))

print(("Wrote %i particles to " % totnevt) + options.output)

outfile.close()
