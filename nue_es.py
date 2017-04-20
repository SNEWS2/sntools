#!/usr/bin/python

from optparse import OptionParser
from math import pi, sin, cos, sqrt, gamma, exp, log
from scipy import integrate, interpolate
from scipy.special import spence
import numpy as np

pid = {"pi0":111, "pi+":211, "k0l":130, "k0s":310, "k+":321,
       "e+":-11, "mu+":-13, "tau+":-15, 
       "nue":12, "nuebar":-12, 
       "numu":14, "numubar":-14, 
       "nutau":16, "nutaubar":-16,
       "p+":2212, "n0":2112}

#holds detector [radius, height] in cm
detectors = {"SuperK":[3368.15/2., 3620., 2.1e+33],
             "HyperK":[7080./2., 5480., 1.4e+34],
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

optdefault = "simData.txt"
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

optchoices = list(detectors.keys())
optdefault = "SuperK"
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

# return direction of an electron with the given energy
#TODO alter for direction of electron rather than ibd proton: approximately the same as the neutrino direction
def direction(energy):
	eneNu = a # TODO find best way to define eneNu
	pMax = 0
	cosT = 0
	nCosTBins = 1000
	cosTBinWidth = 2./nCosTBins
#find the largest value of p
	for j in range(nCosTBins):
		cosT = -1 + cosTBinWidth*(j+0.5) # 1000 steps in the interval [-1,1]
		p = dir_nue_p_sv(eneNu, cosT)
		if p > pMax:
			pMax = p

	while (True):
		cosT = 2*np.random.random() - 1 # randomly distributed in interval [-1,1)
		if dir_nue_p_sv(eneNu, cosT) > pMax*np.random.random():
			sinT = sin(np.arccos(cosT))
			phi = 2 * pi * np.random.random() - pi # randomly distributed in [-pi, pi)
			break

	return (sinT*cos(phi), sinT*sin(phi), cosT)

# probability distribution for the angle at which the positron is emitted
# numerical values are from Ishino-san's code for SK, based on email from Vissani
def dir_nue_p_sv(eneNu, cosT):
	def dir_f1(eneNu):
		return -0.05396 + 0.35824 * (eneNu/100) + 0.03309 * (eneNu/100)**2
	def dir_f2(eneNu):
		return  0.00050 - 0.02390 * (eneNu/100) + 0.14537 * (eneNu/100)**2

	return 0.5 + dir_f1(eneNu) * cosT + dir_f2(eneNu) * (cosT**2 -1./3)

typestr = options.type.replace("+", "plus").replace("-", "minus")

nevtValues=[]
tValues=[]
aValues=[]
xValues=[]
totnevt = 0

#define variables
nE = 9 * detectors[options.detector][2] #number of electrons in detector volume
sin2theta_w = 0.2317 # TODO write in definition 
dSquared = (1.563738e+33)**2
mE = 0.5109989 #MeV
gF=1.16637e-11 #Fermi coupling constant
ro_NC = 1.0126 # +/- 0.0016
eThr = mE # TODO set threshold energy

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
        L=float(L)
        
        #calculate the energy-dependent cross section for (nu_e)-electron scattering        
        
        def l(eNu):
            return sqrt(eNu**2 - mE**2)
        def beta(eNu):
            return l(eNu)/eNu
        def T(eNu, eE):# kinetic energy of recoil electron TODO dependent only on eE
            return eE - mE #total energy of the recoil electron minus the rest mass TODO consider eE_Min must be greater than mE 
        def z(eNu, eE):
            return T(eNu, eE)/eNu # ratio of recoil electron's kinetic energy to incident neutrino energy
        def x(eNu, eE):
            return sqrt(1+(2*mE/T(eNu, eE)))
        def I(eNu, eE):
            return 1/6.0 * (1/3.0 + (3-x(eNu, eE)**2) * ((x(eNu, eE)/2.0)*log((x(eNu, eE)+1)/(x(eNu, eE)-1)) - 1))
        def k(eNu, eE):
            return 0.9791 + 0.0097 * I(eNu, eE) #+/-0.0025 TODO Check whether to include a range to incorporate error
        def gL(eNu, eE):
            return ro_NC * (1/2.0 - k(eNu, eE) * sin2theta_w)-1
        def gR(eNu, eE):
            return -(ro_NC) * k(eNu, eE) * sin2theta_w
           
        #f1 = fMinus(z), f2 = (1-z**2)*fPlus(z), f3 = fMinusPlus(z) for approximate calculation of dSigmadE, see Bahcall, 1995 TODO add web ref

        def f1(eNu, eE):# Non-relativistic approximation
            return ( ((eE/l(eNu)) * log((eE+l(eNu))/mE) - 1)*(2.0*log(1-z(eNu, eE) - mE/(eE+l(eNu))) - log(1.0-z(eNu, eE)) - (1/2.0)*log(z(eNu, eE)) - 5/12.0) 
            + (1/2.0) * (spence(1-beta(eNu)) - spence(1-z(eNu, eE))) - (1/2.0) * (log(1-z(eNu, eE)))**2 
            - (11/12.0 + z(eNu, eE)/2.0) * log(1-z(eNu, eE)) + z(eNu, eE)* (log(z(eNu, eE)) + (1/2.0)*log((2*a)/mE)) 
            - (31/18.0 + (1/12.0)*log(z(eNu, eE)))* beta(eNu) - (11/12.0) * z(eNu, eE) + (z(eNu, eE)**2)/24.0 )
        def f2(eNu, eE):
            return ( ((eE/l(eNu)) * log ((eE + l(eNu))/mE) - 1)* (((1-z(eNu, eE))**2) * (2*log(1 - z(eNu, eE) - (mE/(eE+l(eNu))))-log(1-z(eNu, eE)) - (log (z(eNu, eE)))/2.0 - 2/3.0) 
		 - (z(eNu, eE)**2 * log (z(eNu, eE)) + 1 - z(eNu, eE))/2.0 ) 
            - ((1-z(eNu, eE))**2 / 2.0)*((log(1-z(eNu, eE)))**2 + beta(eNu) * (spence(1-(1-z(eNu, eE))) - log(z(eNu, eE))*log(1-z(eNu, eE)))) 
            + log (1-z(eNu, eE)) * (((z(eNu, eE)**2) / 2.0) * log(z(eNu, eE)) + ((1-z(eNu, eE))/3.0) * (2*z(eNu, eE) - 1/2.0)) 
            - (z(eNu, eE)**2 / 2.0) * spence(1-(1-z(eNu, eE))) - (z(eNu, eE) * (1-2*z(eNu, eE))/3.0) * log (z(eNu, eE)) - z(eNu, eE) * ((1- z(eNu, eE))/6.0) 
            - (beta(eNu)/12.0)* (log(z(eNu, eE)) + (1-z(eNu, eE)) * ((115 - 109 * z(eNu, eE))/6.0)) )
        def f3(eNu, eE):
            return ((eE/l(eNu)) * log ( (eE + l(eNu))/mE) - 1) * 2 * log(1 - z(eNu, eE) - mE/(eE + l(eNu)))
        def dSigmadE(eNu, eE):
            return ( (2*mE*gF**2)/pi * ((gL(eNu, eE)**2) * (1 + (1/137.0/pi) * f1(eNu, eE)) 
            + (gR(eNu, eE)**2) * ((1-z(eNu, eE))**2 + f2(eNu, eE)*(1/137.0/pi))- gR(eNu, eE) * gL(eNu, eE) * (mE/eNu) * z(eNu, eE) * (1+(1/137/pi) * f3(eNu, eE))) )
        
        #calculate the energy-dependent flux per ms        
        alpha = (2*a**2-eNuSquared)/(eNuSquared-a**2)
        def gamma_dist(eNu): #energy distribution of neutrinos
            return (eNu**alpha/gamma(alpha+1))*(((alpha+1)/a)**(alpha+1))*(exp(-(alpha+1)*(eNu/a)))
        def dFluxdE(eNu):
            return 1/(4*pi*dSquared)*((L*624.151)/a)*gamma_dist(eNu)
        
        #calculate the range for electron energy eE TODO calculate the center-of-mass energy of the electron 
        #def s(eNu):
        #    return mE**2 + 2*mE*eNu
        #def eE_cm(eNu):
        #    return sqrt(s(eNu)) # dummy value TODO define centre-of-mass energy

        # Bounds for integration over eE taken from Passera M. (2001) http://arxiv.org/abs/hep-ph/0102212v1
        eE_Min = mE
        def eE_Max(eNu):
            return (mE**2+(2*eNu + mE)**2)/(2*(2*eNu + mE))
        
        #integrate over eE and then eNu to obtain the event rate at time t
        def f(eE, eNu):
            return dSigmadE(eNu, eE)*dFluxdE(eNu)
        def bounds_eNu():
            return [0,50]
        def bounds_eE(eNu):
            return [eE_Min+1, eE_Max(eNu)+1]
        
        #calculate the detector event rate at time t
        simnevt = (nE/0.89) * integrate.nquad(f, [bounds_eE, bounds_eNu]) [0]
        
        #create a list of nevt values at time (t) for input into interpolation function
        nevtValues.append(simnevt)
        
#interpolate the mean energy
interpolatedEnergy = interpolate.pchip(tValues, aValues)

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
    
    if verbose:
    	print ("**************************************")
    	print ("timebin       = %s-%s ms" % (boundsMin, boundsMax))
    	print ("Nevt (theor.) =", binnedNevt)
    	print ("Nevt (actual) =", binnedNevt1ms)
    	print ("mean energy   =", binnedEnergy, "MeV")
    	print ("Now generating events for this bin ...")

    #define particle for each event in time interval
    for i in range(binnedNevt1ms):
        #Define properties of the particle
        t = time - np.random.random()
        ene = np.random.gamma(alpha+1, binnedEnergy/(alpha+1))
        (dirx, diry, dirz) = direction(ene)
        
        #print out [t, pid, energy, dirx, diry, dirz] to file
        outfile.write("%f, -11, %f, %f, %f, %f\n" % (t, ene, dirx, diry, dirz))

print ("**************************************")
print(("Wrote %i particles to " % totnevt) + options.output)

outfile.close()
