# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 12:44:59 2016

@author: liz
"""
from math import pi
import numpy as np
from scipy import integrate

totnevt=0
nP = 4.96*10**34 #number of hydrogen nuclei in detector
d = 10*1.563738*10**32 #distance in MeV
mN = 939.6 #MeV
mP = 938.2 #MeV
mE = 0.511 #MeV
mPi = 139.6 #MeV
delta = mN-mP
mAvg=(mP+mN)/2

#sample mean energy, mean squared energy
eNu, eNuSquared, L = (16, 800, 10**53)
eE = eNu-1.3
sMinusU = (2*mP*(eNu+eE))-(mE**2)
t_eNu_eE = (mN**2) - (mP**2) - (2*mP*(eNu-eE))
f1 = (1-(4.706*(t_eNu_eE/(4*(mAvg**2)))))/((1-(t_eNu_eE/(4*(mAvg**2))))*((1-(t_eNu_eE/710000))**2))
f2 = 3.706/((1-(t_eNu_eE/(4*(mAvg**2))))*((1-(t_eNu_eE/710000))**2))
g1 = (-1.27)/((1-(t_eNu_eE/1000000))**2)
g2 = (2*(mAvg**2)*g1)/((mPi**2) - t_eNu_eE)

A = ((mAvg**2)*((f1**2) - (g1**2))*(t_eNu_eE-(mE**2))) - ((mAvg**2)*(delta**2)*((f1**2)+(g1**2))) - (2*(mE**2)*mAvg*delta*g1*(f1+f2))
B = t_eNu_eE*g1*(f1+f2)
C = ((f1**2) + (g1**2))/4

absMsquared = A-(sMinusU*B)+((sMinusU**2)*C)
dSigmadE = ((((1.16637*(10**(-11)))**2)*(0.9746**2))/(8*pi*(mP**2)*(eNu**2)))*absMsquared*2*mP
a = (eNuSquared-(2*(eNu**2)))/((eNu**2)-eNuSquared)    #for gamma distribution
dFluxdE = (1/(4*pi*(d**2)))*((L/(6.242*(10**5)))/1000*eNu)*np.random.gamma(a+1, eNu/(a+1))#neutrino flux per ms
def f(eE_all, eNu_all):
    return dSigmadE*dFluxdE
def bounds_eNu():
    return [0,100]
def bounds_eE(eNu_all):
    return [0,(eNu_all-1.3)]
simnevt= nP * integrate.nquad(f, [bounds_eE, bounds_eNu])[0] #integrate(dFluxdE*dSigmadE) wrt eE and then eNu
#nevt= np.random.poisson(1, simnevt) #number of events in 1ms interval
#totnevt += len(nevt)

print(dSigmadE, dFluxdE, simnevt)


