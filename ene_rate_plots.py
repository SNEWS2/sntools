#!usr/bin/python
"""
Created on Thu May  4 09:57:20 2017

@author: liz

Plots energy against time and event rate from output of ibd.py or es.py
Useful to check results of event generation phase
"""

import numpy as np
import matplotlib.pyplot as plt

ePartValues=[]
tValues=[]
#import data from ibd.py or es.py output
with open("tmp_es_e.txt") as data: #or tmp_es_e.txt for es
    for line in data:
        t, pid, ePart, xdir, ydir, zdir = line.split(",")
        ePart = float(ePart)
        ePartValues.append(ePart) #create a list of particle energies
        t = float(t)
        tValues.append(t) #create a list of times

simtValues=[]
aValues=[]

#import data from pre-processed supernova simulation data
with open("simData_e.txt") as data2:#or simData_e.txt for es
    for line in data2:
        simt, a, eNuSquared, L = line.split(",")
        simt = float(simt)
        simtValues.append(simt*1000)#create a list of times
        a=float(a)
        aValues.append(a)#create a list of neutrino energies

#create lists of number of events and mean energies for bins of specified widths             
binWidth =2 #time interval in ms
binNr = np.arange(1, 500/binWidth, 1) #time range
meanEnergyValues=[]
nevtValues=[]
   
for i in binNr:
    totalEnergy=0
    nevt=0
    time = 15 + (i*binWidth)
    boundsMin = time - binWidth
    boundsMax = time
    
    for j in range(len(tValues)):
        t = tValues[j]
        ePart = ePartValues[j]
        if boundsMin < t < boundsMax:
            totalEnergy = totalEnergy + ePart
            nevt = nevt+1
    
    if nevt == 0:
        meanEnergyValues.append(totalEnergy)
    else:
        meanEnergy=totalEnergy/nevt
        meanEnergyValues.append(meanEnergy)
    nevtValues.append(nevt) #create a list of number of events in each time bin of specified width

#plot energy against time for charged particles (and neutrinos from supernova simulation)    
plt.title("Charged particle energy") 
plt.ylabel("Energy/MeV")
plt.xlabel("Time/ms")
plt.plot(binWidth*binNr+15, meanEnergyValues)
#plt.plot(simtValues, aValues, 'r')#adds a plot of the neutrino energies direct from the supernova simulation data
plt.show()

#plot event rate
plt.title("Event rate")
plt.ylabel("Events")
plt.xlabel("Time/ms")
plt.plot(binWidth*binNr+15, nevtValues)
plt.show()
