#!usr/bin/python
"""
Created on Thu May  4 09:57:20 2017
Plots energy spectrum of charged particles from neutrino interactions
@author: liz
"""

import matplotlib.pyplot as plt

ePartValues=[]
tValues=[]

#import data from output of es.py/ibd.py
with open("tmp_ibd_eb.txt") as data: #tmp_ibd_eb.txt or tmp_es_e.txt 
    for line in data:
        t, pid, ePart, xdir, ydir, zdir = line.split(",")
        ePart = float(ePart)
        ePartValues.append(ePart)
        t = float(t)
        tValues.append(t)
        
#plot energy spectrum        
plt.hist(ePartValues, bins=50, histtype="step")
plt.title("Energy spectrum")
plt.ylabel("Events")
plt.xlabel("MeV")
plt.show()

#plot energy spectrum with logarithmic y-axis
plt.hist(ePartValues, bins=50, histtype="step", log=(True))
plt.title("Energy spectrum log plot")
plt.ylabel("Events")
plt.xlabel("Energy in MeV")
plt.show()