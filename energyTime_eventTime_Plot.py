#!usr/bin/python
"""
Created on Thu May  4 09:57:20 2017

@author: liz
"""
import numpy as np
import matplotlib.pyplot as plt

ePartValues=[]
tValues=[]

with open("tmp_es_e.txt") as data: #tmp_ibd_eb.txt or tmp_es_e.txt NB must be sorted into ascending order of time
    for line in data:
        t, pid, ePart, xdir, ydir, zdir = line.split(",")
        ePart = float(ePart)
        ePartValues.append(ePart)
        t = float(t)
        tValues.append(t)
        

binWidth =5 #time interval in ms
binNr = np.arange(1, 500/binWidth, 1) #time range
meanEnergyValues=[]
nevtMsValues=[]
   
for i in binNr:
    totalEnergy=0
    nevtMs=0
    time = 15 + (i*binWidth)
    boundsMin = time - binWidth
    boundsMax = time
    
    for j in range(len(tValues)):
        t = tValues[j]
        ePart = ePartValues[j]
        if boundsMin < t < boundsMax:
            totalEnergy = totalEnergy + ePart
            nevtMs = nevtMs+1
    
    if nevtMs == 0:
        meanEnergyValues.append(totalEnergy)
    else:
        meanEnergy=totalEnergy/nevtMs
        meanEnergyValues.append(meanEnergy)
    nevtMsValues.append(nevtMs)
    
#plt.title("") 
plt.ylabel("Mean energy/MeV")
plt.xlabel("Time/ms")
plt.plot(binWidth*binNr+15, meanEnergyValues)
plt.show()

#plt.title("")
plt.ylabel("Event rate")
plt.xlabel("Time/ms")
plt.plot(binWidth*binNr+15, nevtMsValues)
plt.show()

