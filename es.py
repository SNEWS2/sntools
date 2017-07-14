#!/usr/bin/python

from optparse import OptionParser
from math import pi, sin, cos, sqrt, gamma, exp, log
from scipy import integrate, interpolate
#from scipy.special import spence
import numpy as np

pid = {"pi0":111, "pi+":211, "k0l":130, "k0s":310, "k+":321,
       "e+":-11, "mu+":-13, "tau+":-15, 
       "nue":12, "nuebar":-12, 
       "numu":14, "numubar":-14, 
       "nutau":16, "nutaubar":-16,
       "p+":2212, "n0":2112}

#holds detector [radius, height] in cm
detectors = {"SuperK":[3368.15/2., 3620., 2.1e+33],
             "HyperK":[7080./2., 5480., 2*2*1.4e+34],
             "Cylinder_60x74_20inchBandL_14perCent":[7400./2., 6000., 1.7e+34],
             "Cylinder_60x74_20inchBandL_40perCent":[7400./2., 6000., 1.7e+34]}

for pname, no in list(pid.items()):
    if pname.endswith('+'):
        pid[pname.replace('+', '-')] = -1*pid[pname]


parser = OptionParser()

optchoices = list(pid.keys())
optdefault = "e+"
parser.add_option("-t", "--type", dest="type",
                  help="Particle type to be generated. Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
                  metavar="TYPE",
                  choices=optchoices, default=optdefault)

optdefault = "simData_e.txt"
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

optchoices = detectors.keys() #list(detectors.keys()) in python3
optdefault = detectors.keys()[0]
parser.add_option("-d", "--detector", dest="detector",
                  help="Detector configuration. Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
                  choices=optchoices, default=optdefault)

parser.add_option("-v", "--verbose", dest="verbose",
                  help="Verbose output, e.g. for debugging. Off by default.",
                  default=False, action="store_true")

(options, args) = parser.parse_args()

verbose = options.verbose
normalization = float(options.normalization)
if (normalization <= 0 or normalization > 1):
	print("Error: Normalization factor should be in the interval (0,1]. Aborting ...")
	exit()

# return direction of an electron with the given energy in relation to neutrino path
def direction(eneNu):
	
	pMax = 0
	cosT = 0
	nCosTBins = 1000
	cosTBinWidth = 2./nCosTBins
#find the largest value of p
	for j in range(nCosTBins):
		cosT = -1 + cosTBinWidth*(j+0.5) # 1000 steps in the interval [-1,1]
		p = dSigmadCosT(eneNu, cosT)
		if p > pMax:
			pMax = p
	while (True):
		cosT = 2*np.random.random() - 1 # randomly distributed in interval [-1,1)
		if dSigmadCosT(eneNu, cosT) > pMax*np.random.random():
			sinT = sin(np.arccos(cosT))
			phi = 2 * pi * np.random.random() - pi # randomly distributed in [-pi, pi)
			break

	return (sinT*cos(phi), sinT*sin(phi), cosT)

# probability distribution for the scattering angle
# cosT-dependent cross-section can be found at https://www.kvi.nl/~loehner/saf_seminar/2010/neutrino-electron-interactions.pdf
def dSigmadCosT(eneNu, cosT):
	def dir_f1(eneNu, cosT):
		return (sigma0 * 4 * eneNu**2 * (mE+eneNu)**2 * cosT)/((mE+eneNu)**2 - eneNu**2 * cosT**2)**2
	def dir_f2a(eneNu, cosT):
		return 0.73**2 + 0.23**2 * (1 - ((2 * mE * eneNu * cosT**2) / ((mE+eneNu)**2 - eneNu**2 * cosT**2)))**2
	def dir_f2b(eneNu, cosT):
		return (0.73 * 0.23 * 2 * mE**2 * cosT**2)/((mE+eneNu)**2 - (eneNu**2 * cosT**2))
	
	return dir_f1(eneNu, cosT) * (dir_f2a(eneNu, cosT) - dir_f2b(eneNu, cosT)) 

typestr = options.type.replace("+", "plus").replace("-", "minus")

nevtValues=[]
tValues=[]
aValues=[]
eNuSquaredValues=[]
xValues=[]
totnevt = 0
#define variables
nE = 9 * detectors[options.detector][2] #number of electrons in detector volume
sin2theta_w = 0.2317  
dSquared = (1.563738e+33)**2
mE = 0.5109989 #MeV
gF=1.16637e-11 #Fermi coupling constant
ro_NC = 1.0126 # +/- 0.0016
sigma0 = (2*gF**2*mE**2)/pi

#calculate the event rate at each time from the pre-processed data

with open(options.input) as simData:
    if verbose: print("Reading neutrino simulation data from", options.input, "...")
    for line in simData:
        
        #import time, mean energy, mean squared energy and luminosity at time t
        t, a, eNuSquared, L = line.split(",")
        t=float(t)
        t=t*1000
        tValues.append(t)
        a=float(a)
        aValues.append(a)
        eNuSquared = float(eNuSquared)
        eNuSquaredValues.append(eNuSquared)
        L=float(L)
        
        #calculate the energy-dependent cross section for (nu_e)-electron scattering        
        
        def l(eNu):
            return sqrt(eNu**2 - mE**2)
        def beta(eNu):
            return l(eNu)/eNu
        def T(eE):# kinetic energy of recoil electron
            return eE - mE #total energy of the recoil electron minus the rest mass
        def z(eNu, eE):
            return T(eE)/eNu # ratio of recoil electron's kinetic energy to incident neutrino energy
        def x(eE):
            return sqrt(1 + (2*mE/T(eE)))
        def I(eE):
            return 1/6.0 * (1/3.0 + (3 - x(eE)**2) * ((x(eE)/2.0)*log((x(eE) + 1)/(x(eE) - 1)) - 1))
        def k(eE):
            return 0.9791 + 0.0097 * I(eE) #+/-0.0025
        def gL(eE):
            return ro_NC * (1/2.0 - k(eE) * sin2theta_w) - 1
        def gR(eE):
            return -(ro_NC) * k(eE) * sin2theta_w
           
        # f1 = fMinus(z), f2 = (1-z**2)*fPlus(z), f3 = fPlusMinus(z) for approximate calculation of dSigmadE
	# see Bahcall, 1995 DOI:https://doi.org/10.1103/PhysRevD.51.6146 
                   
        def spence(n):
            return integrate.quad(lambda n: log(1-n)/n, 0, n) [0]
        def f1(eNu, eE):
            return ((eE/l(eNu)) * log((eE+l(eNu))/mE) - 1) * (2.0*log(1-z(eNu, eE) - mE/(eE+l(eNu))) - log(1.0-z(eNu, eE)) - (1/2.0)*log(z(eNu, eE)) - 5/12.0) + (1/2.0) * (spence(z(eNu, eE)) - spence(beta(eNu))) - (1/2.0) * (log(1-z(eNu, eE)))**2 - (11/12.0 + z(eNu, eE)/2.0) * log(1-z(eNu, eE)) + z(eNu, eE)* (log(z(eNu, eE)) + (1/2.0)*log((2*a)/mE)) - (31/18.0 + (1/12.0)*log(z(eNu, eE)))* beta(eNu) - (11/12.0) * z(eNu, eE) + (z(eNu, eE)**2)/24.0
        def f2(eNu, eE):
            return ((eE/l(eNu)) * log ((eE + l(eNu))/mE) - 1.) * (((1. - z(eNu, eE))**2) * (2*log(1. - z(eNu, eE) - (mE/(eE+l(eNu))))-log(1.-z(eNu, eE)) - (log (z(eNu, eE)))/2.0 - 2/3.0) - (z(eNu, eE)**2 * log (z(eNu, eE)) + 1 - z(eNu, eE))/2.0 ) - ((1-z(eNu, eE))**2 / 2.0)*((log(1-z(eNu, eE)))**2 + beta(eNu) * (spence(1-z(eNu, eE)) - log(z(eNu, eE))*log(1-z(eNu, eE)))) + log (1-z(eNu, eE)) * (((z(eNu, eE)**2) / 2.0) * log(z(eNu, eE)) + ((1 - z(eNu, eE))/3.0) * (2*z(eNu, eE) - 1/2.0)) - (z(eNu, eE)**2 / 2.0) * spence(1-z(eNu, eE)) - (z(eNu, eE) * (1-2*z(eNu, eE))/3.0) * log (z(eNu, eE)) - z(eNu, eE) * ((1- z(eNu, eE))/6.0) - (beta(eNu)/12.0)* (log(z(eNu, eE)) + (1 - z(eNu, eE)) * ((115 - 109 * z(eNu, eE))/6.0))
        def f3(eNu, eE):
            return ((eE/l(eNu)) * log ( (eE + l(eNu))/mE) - 1) * 2 * log(1 - z(eNu, eE) - mE/(eE + l(eNu)))
        def dSigmadT(eNu, eE):
            return (2*mE*gF**2)/pi * (gL(eE)**2 * (1 + (1/137.0/pi) * f1(eNu, eE)) + gR(eE)**2 * ((1-z(eNu, eE))**2 + f2(eNu, eE)*((1/137.0)/pi))- gR(eE) * gL(eE) * (mE/eNu) * z(eNu, eE) * (1 + ((1/137.)/pi) * f3(eNu, eE)))
            		
        #calculate the energy-dependent flux per ms        
        alpha = (2*a**2-eNuSquared)/(eNuSquared-a**2)
        def gamma_dist(eNu): #energy distribution of neutrinos
            return (eNu**alpha/gamma(alpha + 1))*(((alpha + 1)/a)**(alpha + 1))*(exp(-(alpha + 1)*(eNu/a)))
        def dFluxdE(eNu):
            return 1/(4*pi*dSquared)*((L*624.151)/a)*gamma_dist(eNu)
        
        # Bounds for integration over eE taken from Super-K code
        eE_Min = 0.7 # Cherenkov threshold
        def eE_Max(eNu):
            return ((2*eNu**2)/(2*eNu + mE)) + mE # also explained at //www.kvi.nl/~loehner/saf_seminar/2010/neutrino-electron-interactions.pdf; also eE_Max(cosT) = (2*mE)/(arccos(cosT))**2
        
        #integrate over eE and then eNu to obtain the event rate at time t
        def eNu_min(eE):
            return (T(eE)/2.)*(1 + sqrt(1 + (2*mE)/T(eE))) #from SuperK code, also explained at //www.kvi.nl/~loehner/saf_seminar/2010/neutrino-electron-interactions.pdf
        def f(eE, eNu):
            if eNu>=eNu_min(eE):
            	return dSigmadT(eNu, eE)*dFluxdE(eNu)
            else:
                return 0
        def bounds_eNu():
            return [eNu_min(0.9), 50]#actual minimum 0.33370763820956262 at eE=0.7 MeV calculated using eNu_min(eE) equation
        def bounds_eE(eNu):
            return [eE_Min + 1, eE_Max(eNu) + 1]
        
        #calculate the detector event rate at time t
        simnevt = nE * integrate.nquad(f, [bounds_eE, bounds_eNu]) [0]
        
        #create a list of nevt values at time (t) for input into interpolation function
        nevtValues.append(simnevt)
        
#interpolate the mean energy and mean squared energy
interpolatedEnergy = interpolate.pchip(tValues, aValues)
interpolatedMSEnergy = interpolate.pchip(tValues, eNuSquaredValues)
#interpolate the event rate            
interpolatedNevt = interpolate.pchip(tValues, nevtValues) 

#specify bin width and number of bins for binning to 1ms intervals
binWidth = 1 #time interval in ms
binNr = np.arange(1, 535, 1) #time range

outfile = open(options.output, 'w')
#integrate event rate and energy over each bin
for i in binNr:
    time = 15 + (i*binWidth)
    boundsMin = time - 1
    boundsMax = time

    # calculate expected number of events in this bin and multiply with a factor
    # (1, sin^2(theta_12), cos^2(theta_12)) to take neutrino oscillations into account
    binnedNevt = integrate.quad(interpolatedNevt, boundsMin, boundsMax)[0] * normalization
    #create a poisson distribution of number of events for each bin:
    binnedNevtPoisson = np.random.poisson(binnedNevt, size=1000)
    #randomly select number of events from the Poisson distribution to give the
    #final value for the chosen interval:
    binnedNevt1ms = np.random.choice(binnedNevtPoisson)
    #find the total number of events over all bins
    totnevt += binnedNevt1ms

    binnedEnergy = integrate.quad(interpolatedEnergy, boundsMin, boundsMax)[0]
    binnedMSEnergy = integrate.quad(interpolatedMSEnergy, boundsMin, boundsMax)[0]    
    
    if verbose:
    	print "**************************************"
    	print "timebin       = %s-%s ms" % (boundsMin, boundsMax)
    	print "Nevt (theor.) =", binnedNevt
    	print "Nevt (actual) =", binnedNevt1ms
    	print "mean energy   =", binnedEnergy, "MeV"
    	print "Now generating events for this bin ..."
    
    #define particle for each event in time interval
    for i in range(binnedNevt1ms):
        #Define properties of the particle
        t = time - np.random.random()
        alpha_binned = (2*binnedEnergy**2 - binnedMSEnergy)/(binnedMSEnergy - binnedEnergy**2)
        eneNu = np.random.gamma(alpha_binned + 1, binnedEnergy/(alpha_binned + 1))
        (dirx, diry, dirz) = direction(eneNu)
        while(True): 
            ene = mE + ((2 * mE * eneNu**2 * dirz**2) / ((mE + eneNu)**2 - (eneNu**2 * dirz**2)))
            if ene < (2*eneNu**2/(mE + 2*eneNu))+mE:
                break
        #print out [t, pid, energy, dirx, diry, dirz] to file
        outfile.write("%f, 11, %f, %f, %f, %f\n" % (t, ene, dirx, diry, dirz))

print "**************************************"
print("Wrote %i particles to " % totnevt) + options.output

outfile.close()
