#!/usr/bin/python

from optparse import OptionParser
from math import pi, sin, cos, sqrt, gamma, exp, floor, ceil, log
from scipy import integrate, interpolate
import numpy as np


'''
Setup section.
* Define command line options.
* Parse input options (and perform some sanity checks).
* Read in data from input file.
'''
parser = OptionParser()

optdefault = "infile_e.txt"
parser.add_option("-i", "--input", dest="input",
                  help="Name of the input file. Default: '%s'." \
                      % (optdefault),
                  metavar="FILENAME",
                  default=optdefault)

optdefault = "tmp_es_e.txt"
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

optdefault = "e"
parser.add_option("-f", "--flavor", dest="flavor",
                  help="Neutrino flavor. Choices: 'e', 'eb', 'x', 'xb'. Default: '%s'." \
                      % (optdefault),
                  metavar="flavor",
                  default=optdefault)

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
normalization = float(options.normalization)
if (normalization <= 0 or normalization > 2):
	print("Error: Normalization factor should be in the interval (0,2]. Aborting ...")
	exit()

flavor = options.flavor
if flavor not in ("e", "eb", "x", "xb"):
	print("Error: flavor needs to be one of 'e', 'eb', 'x', 'xb'. Aborting ...")
	exit()

# read data from input file, remove lines with comments and empty lines
with open(options.input) as infile:
	if verbose: print "Reading neutrino simulation data from", options.input, "..."
	indata = [line for line in infile if not (line.startswith("#") or line.isspace())]

# if start time and end time are not given as command line arguments, get them from 1st/last line of input file
_starttime = indata[0].split(",")[0]
_endtime = indata[-1].split(",")[0]
starttime = options.starttime
endtime = options.endtime

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

Based on the Appendices of Bahcall et al. 1995 (https://doi.org/10.1103/PhysRevD.51.6146).
This calculation includes radiative corrections from QCD & QED effects.

For the differences between neutrinos/antineutrinos and a derivation of directionality,
see https://www.kvi.nl/~loehner/saf_seminar/2010/neutrino-electron-interactions.pdf
Careful: This uses different conventions (e.g. minus signs) from Bahcall et al.!
'''
sin2theta_w = 0.2317 # weak mixing angle
mE = 0.5109989 # electron mass (MeV)
gF=1.16637e-11 # Fermi coupling constant
rho_NC = 1.0126 # numerical factor from Bahcall et al.

def dSigmadT(eNu, eE):
	return (2*mE*gF**2)/pi * (gL(eE)**2 * (1 + (1/137.0/pi) * f1(eNu, eE)) + gR(eE)**2 * ((1-z(eNu, eE))**2 + f2(eNu, eE)*((1/137.0)/pi))- gR(eE) * gL(eE) * (mE/eNu) * z(eNu, eE) * (1 + ((1/137.)/pi) * f3(eNu, eE)))

# Appendix A: Radiative Corrections
def l(eE):
	return sqrt(eE**2 - mE**2)
def beta(eNu, eE):
	return l(eE)/eNu
def T(eE):# kinetic energy of recoil electron: total energy minus rest mass
	return eE - mE
def z(eNu, eE):
	return T(eE)/eNu
def x(eE):
	return sqrt(1 + 2*mE/T(eE))
def I(eE):
	return 1/6.0 * (1/3.0 + (3 - x(eE)**2) * (x(eE)/2.0 *log((x(eE) + 1)/(x(eE) - 1)) - 1))
def k(eE):
	if flavor in ("e", "eb"):
		return 0.9791 + 0.0097 * I(eE) #+/-0.0025
	if flavor in ("x", "xb"):
		return 0.9970 - 0.00037 * I(eE) #+/-0.0025

def g1(eE):
	return rho_NC * (1/2.0 - k(eE) * sin2theta_w)
def g2(eE):
	return -rho_NC * k(eE) * sin2theta_w

def gL(eE):
	if   flavor == "e":  return g1(eE) - 1
	elif flavor == "eb": return g2(eE)
	elif flavor == "x":  return g1(eE)
	elif flavor == "xb": return g2(eE)

def gR(eE):
	if   flavor == "e":  return g2(eE)
	elif flavor == "eb": return g1(eE) - 1
	elif flavor == "x":  return g2(eE)
	elif flavor == "xb": return g1(eE)

# Appendix B: QED Effects
# f1 = fMinus(z), f2 = (1-z**2)*fPlus(z), f3 = fPlusMinus(z)
def spence(n):
	return integrate.quad(lambda t: log(abs(1-t))/t, 0, n) [0]
def f1(eNu, eE):
	return ((eE/l(eE)) * log((eE+l(eE))/mE) - 1) * (2.0*log(1-z(eNu, eE) - mE/(eE+l(eE))) - log(1.0-z(eNu, eE)) - (1/2.0)*log(z(eNu, eE)) - 5/12.0) + (1/2.0) * (spence(z(eNu, eE)) - spence(beta(eNu, eE))) - (1/2.0) * (log(1-z(eNu, eE)))**2 - (11/12.0 + z(eNu, eE)/2.0) * log(1-z(eNu, eE)) + z(eNu, eE)* (log(z(eNu, eE)) + (1/2.0)*log((2*eNu)/mE)) - (31/18.0 + (1/12.0)*log(z(eNu, eE)))* beta(eNu, eE) - (11/12.0) * z(eNu, eE) + (z(eNu, eE)**2)/24.0
def f2(eNu, eE):
	return ((eE/l(eE)) * log ((eE + l(eE))/mE) - 1.) * (((1. - z(eNu, eE))**2) * (2*log(1. - z(eNu, eE) - (mE/(eE+l(eE))))-log(1.-z(eNu, eE)) - (log (z(eNu, eE)))/2.0 - 2/3.0) - (z(eNu, eE)**2 * log (z(eNu, eE)) + 1 - z(eNu, eE))/2.0 ) - ((1-z(eNu, eE))**2 / 2.0)*((log(1-z(eNu, eE)))**2 + beta(eNu, eE) * (spence(1-z(eNu, eE)) - log(z(eNu, eE))*log(1-z(eNu, eE)))) + log (1-z(eNu, eE)) * (((z(eNu, eE)**2) / 2.0) * log(z(eNu, eE)) + ((1 - z(eNu, eE))/3.0) * (2*z(eNu, eE) - 1/2.0)) - (z(eNu, eE)**2 / 2.0) * spence(1-z(eNu, eE)) - (z(eNu, eE) * (1-2*z(eNu, eE))/3.0) * log (z(eNu, eE)) - z(eNu, eE) * ((1- z(eNu, eE))/6.0) - (beta(eNu, eE)/12.0)* (log(z(eNu, eE)) + (1 - z(eNu, eE)) * ((115 - 109 * z(eNu, eE))/6.0))
def f3(eNu, eE):
	return ((eE/l(eE)) * log ( (eE + l(eE))/mE) - 1) * 2 * log(1 - z(eNu, eE) - mE/(eE + l(eE)))


# energy of electron scattered into direction cosT by a neutrino with energy eNu
def eneE(eNu, cosT):
	return mE + (2 * mE * eNu**2 * cosT**2) / ((mE + eNu)**2 - eNu**2 * cosT**2)

# distribution of scattering angles
def dSigmadCosT(eNu, cosT):
	dTdCosT = 4 * mE * eNu**2 * (mE+eNu)**2 * cosT / ((mE+eNu)**2 - eNu**2 * cosT**2)**2
	eE = eneE(eNu, cosT)
	return dTdCosT * dSigmadT(eNu, eE)


'''
Astrophysics section.
* neutrinos are well described by a Gamma distribution
* calculate energy-dependent flux at a fiducial distance of 10 kpc
'''
# Convert fiducial distance of 10 kpc into units of MeV**(-1)
# see http://www.wolframalpha.com/input/?i=10+kpc%2F(hbar+*+c)+in+MeV%5E(-1)
dSquared = (1.563738e+33)**2 

# energy distribution of neutrinos
def gamma_dist(eNu, alpha, a):
	return eNu**alpha / gamma(alpha + 1) * ((alpha + 1)/a)**(alpha + 1) * exp(-(alpha + 1)* eNu/a)

def dFluxdE(eNu, luminosity, alpha, a):
	return 1/(4*pi*dSquared) * luminosity/a * gamma_dist(eNu, alpha, a)


'''
Preparation section.
* Parse input data.
* For each time step in the input data, calculate instantaneous event rate.
* Interpolate to get event rate as a function of time.
'''
# double differential event rate
def ddEventRate(eE, eNu, alpha, a):
	if eNu>=eNu_min(eE):
		return dSigmadT(eNu, eE)*dFluxdE(eNu, L, alpha, a)
	else:
		return 0

# Bounds for integration over eE
eE_min = 0.77 # Cherenkov threshold in water (refraction index n=1.34)
def eE_max(eNu):
	return ((2*eNu**2)/(2*eNu + mE)) + mE # this is eneE(eNu, cosT=1); also eE_max(cosT) = (2*mE)/(arccos(cosT))**2
def bounds_eE(eNu, *args): # ignore additional arguments handed over by integrate.nquad()
	return [eE_min, eE_max(eNu)]

# Bounds for integration over eNu
def eNu_min(eE):
	return (T(eE)/2.)*(1 + sqrt(1 + (2*mE)/T(eE))) # inversion of eE_max(eNu)
eNu_max = 50
bounds_eNu = [eNu_min(eE_min), eNu_max]

nevtValues=[]
tValues=[]
aValues=[]
eNuSquaredValues=[]
totnevt = 0
nE = 5 * detectors[options.detector] # number of electrons in detector (8+1+1 per water molecule, i.e. 5 per hydrogen nucleus)

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
	simnevt = nE * integrate.nquad(ddEventRate, [bounds_eE, bounds_eNu], args=(alpha, a)) [0]
	
	# create a list of nevt values at time t for input into interpolation function
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
# Use rejection sampling to get a value from the distribution dist
def rejection_sample(dist, min_val, max_val, n_bins):
	p_max = 0
	bin_width = float(max_val - min_val) / n_bins
	
	for j in range(n_bins):
		val = min_val + bin_width * (j + 0.5)
		p = dist(val)
		if p > p_max:
			p_max = p
	
	while True:
		val = min_val + (max_val - min_val) * np.random.random()
		if p_max * np.random.random() < dist(val):
			break
	
	return val

# return energy of interacting neutrino
def get_eNu(alpha, a):
	dist = lambda _eNu: integrate.quad(ddEventRate, *bounds_eE(_eNu), args=(_eNu, alpha, a))[0]
	eNu = rejection_sample(dist, *bounds_eNu, n_bins=200)
	return eNu

# return direction of scattered electron, if incoming neutrino moves in z direction
def get_direction(eNu):
	dist = lambda x: dSigmadCosT(eNu, x)
	cosT = rejection_sample(dist, 0, 1, 500)
	sinT = sin(np.arccos(cosT))
	phi = 2 * pi * np.random.random() # randomly distributed in [0, 2 pi)
	return (sinT*cos(phi), sinT*sin(phi), cosT)

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
		# Define properties of the particle
		t = boundsMin + np.random.random() * binWidth
		eNu = get_eNu(binnedAlpha, binnedEnergy)
		(dirx, diry, dirz) = get_direction(eNu)
		ene = eneE(eNu, dirz)
		# print out [t, pid, energy, dirx, diry, dirz] to file
		outfile.write("%f, 11, %f, %f, %f, %f\n" % (t, ene, dirx, diry, dirz))

print(("Wrote %i particles to " % totnevt) + options.output)

outfile.close()
