#!/usr/bin/python

from optparse import OptionParser
from math import pi, sin, cos, sqrt, gamma, exp, floor, ceil, log
from scipy import integrate, interpolate
import numpy as np

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


# return direction of an electron with the given energy in relation to neutrino path
def direction(eNu):
	pMax = 0
	nCosTBins = 500
	cosTBinWidth = 1./nCosTBins
	# find the largest value of p
	for j in range(nCosTBins):
		cosT = cosTBinWidth * (j+0.5) # 500 steps in the interval [0,1)
		p = dSigmadCosT(eNu, cosT)
		if p > pMax:
			pMax = p
	while (True):
		cosT = 2*np.random.random() - 1 # randomly distributed in interval [-1,1)
		if dSigmadCosT(eNu, cosT) > pMax * np.random.random():
			sinT = sin(np.arccos(cosT))
			phi = 2 * pi * np.random.random() # randomly distributed in [0, 2 pi)
			break
	return (sinT*cos(phi), sinT*sin(phi), cosT)

def eneE(eNu, cosT):
	return mE + (2 * mE * eNu**2 * cosT**2) / ((mE + eNu)**2 - eNu**2 * cosT**2)

# probability distribution for the scattering angle, see derivation at
# https://www.kvi.nl/~loehner/saf_seminar/2010/neutrino-electron-interactions.pdf
def dSigmadCosT(eNu, cosT):
	dTdCosT = 4 * mE * eNu**2 * (mE+eNu)**2 * cosT / ((mE+eNu)**2 - eNu**2 * cosT**2)**2
	eE = eneE(eNu, cosT)
	return dTdCosT * dSigmadT(eNu, eE)

nevtValues=[]
tValues=[]
aValues=[]
eNuSquaredValues=[]
xValues=[]
totnevt = 0
# define variables
nE = 5 * detectors[options.detector] # number of electrons in detector volume (8+1+1 per water molecule, i.e. 5 per hydrogen nucleus)
sin2theta_w = 0.2317
dSquared = (1.563738e+33)**2
mE = 0.5109989 # MeV
gF=1.16637e-11 # Fermi coupling constant
rho_NC = 1.0126 # +/- 0.0016

# calculate the event rate at each time from the pre-processed data
for line in indata:
	# import time, mean energy, mean squared energy and luminosity at time t
	t, a, eNuSquared, L = line.split(",")
	t=float(t) * 1000 # convert to ms
	tValues.append(t)
	a=float(a)
	aValues.append(a)
	eNuSquared = float(eNuSquared)
	eNuSquaredValues.append(eNuSquared)
	L=float(L)
	
	# calculate the energy-dependent cross section for (nu_e)-electron scattering
	# see Bahcall, 1995 DOI:https://doi.org/10.1103/PhysRevD.51.6146 (Appendices)
	
	# Appendix A: Radiative Corrections
	def l(eNu):
		return sqrt(eNu**2 - mE**2)
	def beta(eNu):
		return l(eNu)/eNu
	def T(eE):# kinetic energy of recoil electron
		return eE - mE # total energy of the recoil electron minus the rest mass
	def z(eNu, eE):
		return T(eE)/eNu # ratio of recoil electron's kinetic energy to incident neutrino energy
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
		return integrate.quad(lambda n: log(1-n)/n, 0, n) [0]
	def f1(eNu, eE):
		return ((eE/l(eNu)) * log((eE+l(eNu))/mE) - 1) * (2.0*log(1-z(eNu, eE) - mE/(eE+l(eNu))) - log(1.0-z(eNu, eE)) - (1/2.0)*log(z(eNu, eE)) - 5/12.0) + (1/2.0) * (spence(z(eNu, eE)) - spence(beta(eNu))) - (1/2.0) * (log(1-z(eNu, eE)))**2 - (11/12.0 + z(eNu, eE)/2.0) * log(1-z(eNu, eE)) + z(eNu, eE)* (log(z(eNu, eE)) + (1/2.0)*log((2*a)/mE)) - (31/18.0 + (1/12.0)*log(z(eNu, eE)))* beta(eNu) - (11/12.0) * z(eNu, eE) + (z(eNu, eE)**2)/24.0
	def f2(eNu, eE):
		return ((eE/l(eNu)) * log ((eE + l(eNu))/mE) - 1.) * (((1. - z(eNu, eE))**2) * (2*log(1. - z(eNu, eE) - (mE/(eE+l(eNu))))-log(1.-z(eNu, eE)) - (log (z(eNu, eE)))/2.0 - 2/3.0) - (z(eNu, eE)**2 * log (z(eNu, eE)) + 1 - z(eNu, eE))/2.0 ) - ((1-z(eNu, eE))**2 / 2.0)*((log(1-z(eNu, eE)))**2 + beta(eNu) * (spence(1-z(eNu, eE)) - log(z(eNu, eE))*log(1-z(eNu, eE)))) + log (1-z(eNu, eE)) * (((z(eNu, eE)**2) / 2.0) * log(z(eNu, eE)) + ((1 - z(eNu, eE))/3.0) * (2*z(eNu, eE) - 1/2.0)) - (z(eNu, eE)**2 / 2.0) * spence(1-z(eNu, eE)) - (z(eNu, eE) * (1-2*z(eNu, eE))/3.0) * log (z(eNu, eE)) - z(eNu, eE) * ((1- z(eNu, eE))/6.0) - (beta(eNu)/12.0)* (log(z(eNu, eE)) + (1 - z(eNu, eE)) * ((115 - 109 * z(eNu, eE))/6.0))
	def f3(eNu, eE):
		return ((eE/l(eNu)) * log ( (eE + l(eNu))/mE) - 1) * 2 * log(1 - z(eNu, eE) - mE/(eE + l(eNu)))
	
	# Complete formula
	def dSigmadT(eNu, eE):
		return (2*mE*gF**2)/pi * (gL(eE)**2 * (1 + (1/137.0/pi) * f1(eNu, eE)) + gR(eE)**2 * ((1-z(eNu, eE))**2 + f2(eNu, eE)*((1/137.0)/pi))- gR(eE) * gL(eE) * (mE/eNu) * z(eNu, eE) * (1 + ((1/137.)/pi) * f3(eNu, eE)))
	
	# calculate the energy-dependent flux per ms
	alpha = (2*a**2-eNuSquared)/(eNuSquared-a**2)
	def gamma_dist(eNu): # energy distribution of neutrinos
		return (eNu**alpha/gamma(alpha + 1))*(((alpha + 1)/a)**(alpha + 1))*(exp(-(alpha + 1)*(eNu/a)))
	def dFluxdE(eNu):
		return 1/(4*pi*dSquared)*((L*624.151)/a)*gamma_dist(eNu)
	
	# Bounds for integration over eE taken from Super-K code
	eE_Min = 0.7 # Cherenkov threshold
	def eE_Max(eNu):
		return ((2*eNu**2)/(2*eNu + mE)) + mE # also explained at //www.kvi.nl/~loehner/saf_seminar/2010/neutrino-electron-interactions.pdf; also eE_Max(cosT) = (2*mE)/(arccos(cosT))**2
	
	# integrate over eE and then eNu to obtain the event rate at time t
	def eNu_min(eE):
		return (T(eE)/2.)*(1 + sqrt(1 + (2*mE)/T(eE))) # from SuperK code, also explained at //www.kvi.nl/~loehner/saf_seminar/2010/neutrino-electron-interactions.pdf
	def f(eE, eNu):
		if eNu>=eNu_min(eE):
			return dSigmadT(eNu, eE)*dFluxdE(eNu)
		else:
			return 0
	def bounds_eNu():
		return [eNu_min(0.9), 50]# actual minimum 0.33370763820956262 at eE=0.7 MeV calculated using eNu_min(eE) equation
	def bounds_eE(eNu):
		return [eE_Min + 1, eE_Max(eNu) + 1]
	
	# calculate the detector event rate at time t
	simnevt = nE * integrate.nquad(f, [bounds_eE, bounds_eNu]) [0]
	
	# create a list of nevt values at time t for input into interpolation function
	nevtValues.append(simnevt)

# interpolate the mean energy and mean squared energy
interpolatedEnergy = interpolate.pchip(tValues, aValues)
interpolatedMSEnergy = interpolate.pchip(tValues, eNuSquaredValues)

# interpolate the event rate
interpolatedNevt = interpolate.pchip(tValues, nevtValues)

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
		eNu = np.random.gamma(binnedAlpha + 1, binnedEnergy/(binnedAlpha + 1))
		(dirx, diry, dirz) = direction(eNu)
		ene = eneE(eNu, dirz)
		# print out [t, pid, energy, dirx, diry, dirz] to file
		outfile.write("%f, 11, %f, %f, %f, %f\n" % (t, ene, dirx, diry, dirz))

print(("Wrote %i particles to " % totnevt) + options.output)

outfile.close()
