#!/usr/bin/python

from optparse import OptionParser
import random
from math import pi, sin, cos, sqrt, gamma, exp
from scipy import integrate, interpolate
import numpy as np

pid = {"pi0":111, "pi+":211, "k0l":130, "k0s":310, "k+":321,
       "e+":-11, "mu+":-13, "tau+":-15, 
       "nue":12, "nuebar":-12, 
       "numu":14, "numubar":-14, 
       "nutau":16, "nutaubar":-16,
       "p+":2212, "n0":2112}

#holds detector [radius, height] in cm
detectors = {"SuperK":[3368.15/2., 3620.],
             "Cylinder_60x74_20inchBandL_14perCent":[7400./2., 6000.],
             "Cylinder_60x74_20inchBandL_40perCent":[7400./2., 6000.]}

for pname, no in list(pid.items()):
    if pname.endswith('+'):
        pid[pname.replace('+', '-')] = -1*pid[pname]





parser = OptionParser()
optdefault = 1
parser.add_option("-N", "--nfiles", dest="nfiles",
                  help="number of files of particles to produce. Default: %s" \
                      % (optdefault),
                  metavar="#", default=optdefault)
optchoices = list(pid.keys())
optdefault = "e+"
parser.add_option("-t", "--type", dest="type",
                  help="Particle type to be generated. Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
                  metavar="TYPE",
                  choices=optchoices, default=optdefault)
optdefault = 1000
parser.add_option("-e", "--energy", dest="energy",
                  help="Particle energy to be generated in MeV. Default: %s" \
                      % (optdefault),
                  metavar="ENERGY", default=optdefault)
optchoices = ["center", "random", "minusx", "plusx", "minusz", "plusz"]
optdefault = optchoices[1]
parser.add_option("-v", "--vertex", dest="vertname",
                  help="Type of vertex. Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
                  choices=optchoices, default=optdefault)
optchoices = ["4pi", "towall", "tocap"]
optdefault = optchoices[0]
parser.add_option("-d", "--direction", dest="dirname",
                  help="Type of direction. Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
                  choices=optchoices, default=optdefault)
optchoices = list(detectors.keys())
optdefault = "SuperK"
parser.add_option("-w", "--detector", dest="detector",
                  help="Detector water volume to use (for vertex positioning). Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
                  choices=optchoices, default=optdefault)

(options, args) = parser.parse_args()

options.vertname = options.vertname.lower()
options.dirname = options.dirname.lower()


nfiles = int(options.nfiles)

def partPrint(p, f, recno):
    f.write("$ begin\n")
    f.write("$ nuance 0\n")
    rad    = detectors[options.detector][0] - 20.
    height = detectors[options.detector][1] - 20.
    while True:
        x = random.uniform(-rad, rad)
        y = random.uniform(-rad, rad)
        if x**2 + y**2 < rad**2: break
    z = random.uniform(-height/2, height/2)
    f.write("$ vertex %.5f %.5f %.5f %.5f\n" % (x, y, z, p["time"]))
    printTrack(nu, f, -1)   # "Neutrino" Track
    printTrack(prot, f, -1) # "Target" track
    f.write("$ info 0 0 %i\n" % recno)
    th = random.random()*2*pi
    u = 1.-2*random.random()
    x = sqrt(1.-u**2)*cos(th)
    y = sqrt(1.-u**2)*sin(th)
    z = u
    p["direction"] = (x, y, z)
    #th = random.random()*pi
    #phi = random.random()*2*pi
    #p["direction"] = (cos(phi)*cos(th), sin(phi)*cos(th), sin(th))
       
    printTrack(p, f)    # Outgoing Particle Track
    f.write("$ end\n")

def printTrack(p, f, code=0):
    f.write("$ track %(type)i %(energy).5f " % p)
    f.write("%.5f %.5f %.5f %i\n" % (p["direction"]+(code,)))

for fileno in range(nfiles):
            typestr = options.type.replace("+", "plus").replace("-", "minus")
            
            filename="%s_%s_%s_%s_%03i.kin" % (typestr, options.vertname, options.dirname, options.detector, fileno)
        
            outfile = open(filename, 'w')
nevtValues=[]
tValues=[]
aValues=[]
totnevt = 0
#define variables
nP_HK = 4.96e+34 #number of hydrogen nuclei in whole detector volume; needs to be updated for design changes
dSquared = (1.563738e+33)**2
mN = 939.5654 #MeV
mP = 938.2721 #MeV
mE = 0.5109989 #MeV
mPi = 139.57018 #MeV
delta = mN-mP
mAvg=(mP+mN)/2
gF=1.16637e-11 #Fermi coupling constant
eThr=((mN+mE)**2 - mP**2)/(2*mP) #threshold energy for IBD

#calculate the event rate at each time from the pre-processed data

with open('simData.txt') as simData:
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

        #calculate the energy-dependent cross section for IBD        
        
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
           
        #AM, BM and CM for approximate calculation of absMsquared, AM1, BM1,
        #CM1 for more precise calculation
        def AM(eNu, eE):
            return (mAvg**2 * (f1(eNu, eE)**2 - g1(eNu, eE)**2) *
            (t_eNu_eE(eNu, eE)-mE**2)) - (mAvg**2 * delta**2 * 
            (f1(eNu, eE)**2 + g1(eNu, eE)**2)) - (2 * mE**2 * mAvg * delta * g1(eNu, eE) *
            (f1(eNu, eE)+f2(eNu, eE)))
        def AM1(eNu, eE):
            return  ((1/16)*(t_eNu_eE(eNu, eE)-mE**2) * ((4*(f1(eNu, eE)**2)*
            ((4*mAvg**2)+t_eNu_eE(eNu, eE)+mE**2))
            +(4*(g1(eNu, eE)**2)*(-4*(mAvg**2)+t_eNu_eE(eNu, eE)+mE**2))
            +((f2(eNu, eE)**2)*(((t_eNu_eE(eNu, eE)**2)/(mAvg**2))+4*t_eNu_eE(eNu, eE)+4*mE**2))
            +(4*(mE**2)*t_eNu_eE(eNu, eE)*(g2(eNu, eE)**2)/(mAvg**2)) 
            +(8*f1(eNu, eE)*f2(eNu, eE)*(2*t_eNu_eE(eNu, eE)+mE**2))
            +(16*(mE**2)*g1(eNu, eE)*g2(eNu, eE)))
            -(delta**2)*((4*(f1(eNu, eE)**2)+t_eNu_eE(eNu, eE)*(f2(eNu, eE)**2)/mAvg**2)*
            (4*(mAvg**2)+t_eNu_eE(eNu, eE)-mE**2)
            +(4*(g1(eNu, eE)**2)*(4*(mAvg**2)-t_eNu_eE(eNu, eE)+mE**2))
            +4*(mE**2)*(g2(eNu, eE)**2)*(t_eNu_eE(eNu, eE)-mE**2)/(mAvg**2)
            +8*f1(eNu, eE)*f2(eNu, eE)*(2*t_eNu_eE(eNu, eE)-mE**2)
            +16*(mE**2)*g1(eNu, eE)*g2(eNu, eE))
            -32*(mE**2)*mAvg*delta*g1(eNu, eE)*(f1(eNu, eE)+f2(eNu, eE)))
        def BM(eNu, eE):
            return t_eNu_eE(eNu, eE)*g1(eNu, eE)*(f1(eNu, eE)+f2(eNu, eE))
        def BM1(eNu, eE):
            return ((1/16)*(16*t_eNu_eE(eNu, eE)*g1(eNu, eE)*(f1(eNu, eE)+f2(eNu, eE))
            +4*(mE**2)*delta*((f2(eNu, eE)**2)+f1(eNu, eE)*f2(eNu, eE)+2*g1(eNu, eE)*g2(eNu, eE))/mAvg))
        def CM(eNu, eE):
            return ((f1(eNu, eE)**2) + (g1(eNu, eE)**2))/4
        def CM1(eNu, eE):
            return (1/16)*(4*((f1(eNu, eE)**2)+(g1(eNu, eE)**2)) - ((t_eNu_eE(eNu, eE)*(f2(eNu, eE)**2))/(mAvg**2)))
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
        simnevt= (nP_HK/0.89) * integrate.nquad(f, [bounds_eE, bounds_eNu]) [0]
        
        #create a list of nevt values at time (t) for input into interpolation function
        nevtValues.append(simnevt)

#interpolate the mean energy
interpolatedEnergy = interpolate.pchip(tValues, aValues)

#interpolate the event rate            
interpolatedNevt = interpolate.pchip(tValues, nevtValues) 

#specify bin width and number of bins for binning to 1ms intervals
binWidth=1#time interval in ms
binNr = np.arange(1, 535, 1)#time range

#integrate event rate and energy over each bin
for i in binNr:
    
    time=15+(i*binWidth)
    interpolatedNevtValues = interpolatedNevt(time)
    
    interpolatedEnergyValues = interpolatedEnergy(time)    
    
    boundsMin = time-1
    boundsMax = time
    binnedNevt = integrate.quad(interpolatedNevt, boundsMin, boundsMax)[0]
    #create a poisson distribution of number of events for each bin:
    binnedNevtPoisson = np.random.poisson(binnedNevt, size=1000)
    #randomly select number of events from the Poisson distribution to give the
    #final value for the chosen interval:
    binnedNevt1ms = np.random.choice(binnedNevtPoisson)
    #find the total number of events over all bins
    totnevt += binnedNevt1ms

    binnedEnergy = integrate.quad(interpolatedEnergy, boundsMin, boundsMax)[0]
    

    #define particle for each event in time interval
    for i in range(binnedNevt1ms):
        #Define the particle
        particle = {"vertex":(),
                    "time": t, #To do: change this to new values for t
                    "type":pid[options.type],
                    "energy":np.random.gamma(alpha+1, binnedEnergy/(alpha+1)),
                    "direction":()}


        nu =   {"type":pid["numu"], "energy":1000.0, #removed energy+
               "direction":(1, 0, 0)}
        prot = {"type":pid["p+"], "energy":935.9840,
               "direction":(0, 0, 1)}
      
        partPrint(particle, outfile, i)

print(("Writing %i particles to " % totnevt) + filename)

outfile.write("$ stop")
outfile.close()

