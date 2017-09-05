#!/usr/bin/python

from optparse import OptionParser
import random
from math import pi, sin, cos, sqrt, gamma, exp
from scipy import integrate, interpolate
import numpy as np

parser = OptionParser()

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

# number of free protons (i.e. H nuclei) in each detector
detectors = {"SuperK": 2.1e+33,
             "HyperK": 1.4e+34}
optchoices = detectors.keys()
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

# return direction of a positron with the given energy
def direction(energy):
	eneNu = energy + eThr
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

nevtValues=[]
tValues=[]
aValues=[]
totnevt = 0
#define variables
nP = detectors[options.detector] # number of protons in detector volume
dSquared = (1.563738e+33)**2 # 10 kpc in units of MeV**(-1), see http://www.wolframalpha.com/input/?i=10+kpc%2F(hbar+*+c)+in+MeV%5E(-1)
mN = 939.5654 #MeV
mP = 938.2721 #MeV
mE = 0.5109989 #MeV
mPi = 139.57018 #MeV
delta = mN-mP
mAvg=(mP+mN)/2
gF=1.16637e-11 #Fermi coupling constant
eThr=((mN+mE)**2 - mP**2)/(2*mP) #threshold energy for IBD

#calculate the event rate at each time from the pre-processed data

with open(options.input) as simData:
    if verbose: print "Reading neutrino simulation data from", options.input, "..."
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

        # the following code implements the energy-dependent cross section for IBD
        # from Strumia/Vissani (2003), arXiv:astro-ph/0302055
        
        def t_eNu_eE(eNu, eE):
            return mN**2 - mP**2 - 2*mP*(eNu-eE)
        def x(eNu, eE):
            return t_eNu_eE(eNu, eE)/(4*mAvg**2)
        def y (eNu, eE):
            return 1-(t_eNu_eE(eNu, eE)/710000)
        def z (eNu, eE):
            return 1-(t_eNu_eE(eNu, eE)/1000000)
        def f1 (eNu, eE):
            return (1-(4.706*x(eNu, eE)))/((1-x(eNu, eE))*y(eNu, eE)**2)
        def f2 (eNu, eE):
            return 3.706/((1-x(eNu, eE))*y(eNu, eE)**2)
        def g1(eNu, eE):
            return (-1.27)/z(eNu, eE)**2
        def g2(eNu, eE):
            return (2 * g1(eNu, eE) * mAvg**2)/(mPi**2 - t_eNu_eE(eNu, eE)) 
           
        # AM, BM and CM for approximate calculation of absMsquared,
        # AM1, BM1, CM1 for more precise calculation
        def AM(eNu, eE):
            return (mAvg**2 * (f1(eNu, eE)**2 - g1(eNu, eE)**2) *
            (t_eNu_eE(eNu, eE)-mE**2)) - (mAvg**2 * delta**2 * 
            (f1(eNu, eE)**2 + g1(eNu, eE)**2)) - (2 * mE**2 * mAvg * delta * g1(eNu, eE) *
            (f1(eNu, eE)+f2(eNu, eE)))
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
            return (2*mP*(eNu+eE))-mE**2
        def absMsquared(eNu, eE):
            return AM1(eNu, eE)-(sMinusU(eNu, eE)*BM1(eNu, eE))+((sMinusU(eNu, eE)**2)*CM1(eNu, eE))
        def dSigmadE(eNu, eE):
            return (2*mP*gF**2 * (0.9746**2))/(8 * pi * mP**2 * eNu**2)*absMsquared(eNu, eE)
        
        #calculate the energy-dependent flux per ms        
        alpha = (2*a**2-eNuSquared)/(eNuSquared-a**2)
        def gamma_dist(eNu): #energy distribution of neutrinos
            return (eNu**alpha/gamma(alpha+1))*(((alpha+1)/a)**(alpha+1))*(exp(-(alpha+1)*(eNu/a)))
        def dFluxdE(eNu):
            return 1/(4*pi*dSquared)*((L*624.151)/a)*gamma_dist(eNu)
        
        #calculate range for eE from eNu
        def s(eNu):
            return 2*mP*eNu + mP**2
        def pE_cm(eNu):
            return (sqrt((s(eNu)-(mN-mE)**2)*(s(eNu)-(mN+mE)**2)))/(2*sqrt(s(eNu)))
        #def eNu_cm(eNu):
         #   return (s(eNu)-mP**2)/(2*sqrt(s(eNu)))
        def eE_cm(eNu):
            return (s(eNu)-mN**2+mE**2)/(2*sqrt(s(eNu)))
        delta_cm = (mN**2-mP**2-mE**2)/(2*mP)
        def eE_Min(eNu):
            return eNu - delta_cm - (eNu/sqrt(s(eNu)) * (eE_cm(eNu) + pE_cm(eNu)))
        def eE_Max(eNu):
            return eNu - delta_cm - (eNu/sqrt(s(eNu)) *(eE_cm(eNu) - pE_cm(eNu)))
        
        #integrate over eE and then eNu to obtain the event rate at time t
        def f(eE, eNu):
            return dSigmadE(eNu, eE)*dFluxdE(eNu)
        def bounds_eNu():
            return [eThr,50]
        def bounds_eE(eNu):
            return [eE_Min(eNu)+1, eE_Max(eNu)+1]
        
        #calculate the detector event rate at time t
        simnevt = nP * integrate.nquad(f, [bounds_eE, bounds_eNu]) [0]
        
        #create a list of nevt values at time (t) for input into interpolation function
        nevtValues.append(simnevt)

#interpolate the mean energy
interpolatedEnergy = interpolate.pchip(tValues, aValues)

#interpolate the event rate            
interpolatedNevt = interpolate.pchip(tValues, nevtValues) 

#specify bin width and number of bins for binning to 1ms intervals
binWidth = 1 #time interval in ms
binNr = np.arange(1, 535/binWidth, 1) #time range

outfile = open(options.output, 'w')
#integrate event rate and energy over each bin
for i in binNr:
    time = 15 + (i*binWidth)
    boundsMin = time - binWidth
    boundsMax = time

    # calculate expected number of events in this bin and multiply with a factor
    # (1, sin^2(theta_12), cos^2(theta_12)) to take neutrino oscillations into account
    binnedNevt = integrate.quad(interpolatedNevt, boundsMin, boundsMax)[0] * normalization
    # randomly select number of events in this bin from Poisson distribution around binnedNevt:
    binnedNevtRnd = np.random.choice(np.random.poisson(binnedNevt, size=1000))
    #find the total number of events over all bins
    totnevt += binnedNevtRnd

    binnedEnergy = integrate.quad(interpolatedEnergy, boundsMin, boundsMax)[0]
    
    if verbose:
    	print "**************************************"
    	print "timebin       = %s-%s ms" % (boundsMin, boundsMax)
    	print "Nevt (theor.) =", binnedNevt
    	print "Nevt (actual) =", binnedNevtRnd
    	print "mean energy   =", binnedEnergy, "MeV"
    	print "Now generating events for this bin ..."

    #define particle for each event in time interval
    for i in range(binnedNevtRnd):
        #Define properties of the particle
        t = time - np.random.random()
        ene = np.random.gamma(alpha+1, binnedEnergy/(alpha+1))
        (dirx, diry, dirz) = direction(ene)
        
        # print out [t, pid, energy, dirx, diry, dirz] to file
        outfile.write("%f, -11, %f, %f, %f, %f\n" % (t, ene, dirx, diry, dirz))

print "**************************************"
print(("Wrote %i particles to " % totnevt) + options.output)

outfile.close()
