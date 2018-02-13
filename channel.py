#!/usr/bin/python

from importlib import import_module
from math import pi, sin, cos, acos, gamma, exp, floor, ceil
import numpy as np
import random
from scipy import integrate, interpolate

def main(channel="ibd", input="infile_eb.txt", output="tmp_ibd_eb.txt",
         normalization=1.0, detector="SuperK", starttime=None, endtime=None,
         verbose=False):

    '''
    Setup section.
    * Read in data from input file and apply start & end time.
    * Import channel-specific functions & parameters from separate files
    '''
    # inner detector mass, in metric kt
    detectors = {"SuperK": 32.5,
                 "HyperK": 220} # one-tank configuration

    # read data from input file, remove lines with comments and empty lines
    with open(input) as infile:
        if verbose: print "Reading neutrino simulation data from", input, "..."
        raw_indata = [map(float, line.split(",")) for line in infile if not (line.startswith("#") or line.isspace())]

    # if start time and end time are not given as command line arguments, get them from 1st/last line of input file
    _starttime = raw_indata[0][0]
    _endtime = raw_indata[-1][0]

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

    # remove entries outside the requested range
    indata = []
    for (i, entry) in enumerate(raw_indata):
        if i == 0: continue
        if entry[0] > starttime:
            indata.append(raw_indata[i-1])
            if entry[0] > endtime:
                indata.append(entry)
                break

    starttime = ceil(float(starttime) * 1000) # convert to ms
    endtime = floor(float(endtime) * 1000)
    duration = endtime - starttime

    channel_module = import_module("interaction-channels." + channel)
    dSigma_dE = getattr(channel_module,"dSigma_dE")
    get_eE = getattr(channel_module,"get_eE")
    dSigma_dCosT = getattr(channel_module,"dSigma_dCosT")
    bounds_eE = getattr(channel_module,"bounds_eE")
    bounds_eNu = getattr(channel_module,"bounds_eNu")
    targets_per_molecule = getattr(channel_module,"targets_per_molecule")
    pid = getattr(channel_module,"pid")


    '''
    Astrophysics section.
    * neutrinos are well described by a Gamma distribution (arXiv:1211.3920)
    * calculate energy-dependent flux at a fiducial distance of 10 kpc
    '''
    # Convert fiducial distance of 10 kpc into units of MeV**(-1)
    # see http://www.wolframalpha.com/input/?i=10+kpc%2F(hbar+*+c)+in+MeV%5E(-1)
    dSquared = (1.563738e+33)**2

    # energy distribution of neutrinos
    def gamma_dist(eNu, alpha, a):
        return eNu**alpha / gamma(alpha + 1) * ((alpha + 1)/a)**(alpha + 1) * exp(-(alpha + 1)* eNu/a)

    def dFluxdE(eNu, luminosity, alpha, a):
        # The `normalization` factor takes into account the oscillation probability
        # as well as the distance (if not equal to 10 kpc).
        return 1/(4*pi*dSquared) * luminosity/a * gamma_dist(eNu, alpha, a) * normalization


    '''
    Preparation section.
    * Parse input data.
    * For each time step in the input data, calculate instantaneous event rate.
    * Interpolate to get event rate as a function of time.
    '''
    # double differential event rate
    def ddEventRate(eE, eNu, alpha, a, L):
        return dSigma_dE(eNu, eE)*dFluxdE(eNu, L, alpha, a)


    nevtValues = []
    tValues = []
    aValues = []
    eNuSquaredValues = []
    totnevt = 0
    molecules_per_kt = 3.343e+31 # number of water molecules in one kt (assuming 18 g/mol)
    n_targets = targets_per_molecule * molecules_per_kt * detectors[detector]

    for (t, a, eNuSquared, L) in indata:
        # get time, mean energy, mean squared energy, luminosity
        tValues.append(1000 * t) # convert to ms
        aValues.append(a)
        eNuSquaredValues.append(eNuSquared)
        L = L * 624.151 # convert from erg/s to MeV/ms

        alpha = (2*a**2 - eNuSquared) / (eNuSquared - a**2)

        # integrate over eE and then eNu to obtain the event rate at time t
        simnevt = n_targets * integrate.nquad(ddEventRate, [bounds_eE, bounds_eNu], args=(alpha, a, L)) [0]

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
    def rejection_sample(dist, min_val, max_val, n_bins=100):
        p_max = 0
        j_max = 0
        bin_width = float(max_val - min_val) / n_bins

        # Iterative approach to speed up finding the maximum of `dist`.
        # Assumes that `dist` does not oscillate very quickly.
        # First, use coarse binning to find the approximate maximum:
        for j in range(0, n_bins, 10):
            val = min_val + bin_width * (j + 0.5)
            p = dist(val)
            if p > p_max:
                p_max = p
                j_max = j
        # Then, use finer binning around the approximate maximum.
        for j in range(max(j_max-9, 0), min(j_max+10, n_bins)):
            val = min_val + bin_width * (j + 0.5)
            p = dist(val)
            if p > p_max:
                p_max = p

        while True:
            val = min_val + (max_val - min_val) * random.random()
            if p_max * random.random() < dist(val):
                break

        return val

    # return energy of interacting neutrino
    def get_eNu(alpha, a):
        # L is a constant factor so, for the purpose of rejection sampling, we can set it to 1.
        dist = lambda _eNu: integrate.quad(ddEventRate, *bounds_eE(_eNu), args=(_eNu, alpha, a, 1.))[0]
        eNu = rejection_sample(dist, *bounds_eNu, n_bins=200)
        return eNu

    # return direction of scattered electron, if incoming neutrino moves in z direction
    def get_direction(eNu):
        dist = lambda x: dSigma_dCosT(eNu, x)
        cosT = rejection_sample(dist, -1, 1, 200)
        sinT = sin(acos(cosT))
        phi = 2 * pi * random.random() # randomly distributed in [0, 2 pi)
        return (sinT*cos(phi), sinT*sin(phi), cosT)

    binWidth = 1 # bin width in ms
    binNr = np.arange(1, floor(duration/binWidth)+1) # number of full-width bins
    if verbose:
        print "Now generating events in %s ms bins between %s-%s ms" % (binWidth, starttime, endtime)
        print "**************************************"

    outfile = open(output, 'w')
    # integrate event rate and energy over each bin
    for i in binNr:
        boundsMin = starttime + (i-1)*binWidth
        boundsMax = starttime + i*binWidth

        # calculate expected number of events in this bin
        binnedNevt = interpolatedNevt(boundsMin + 0.5*binWidth)
        # randomly select number of events in this bin from Poisson distribution around binnedNevt:
        binnedNevtRnd = np.random.poisson(binnedNevt)
        # find the total number of events over all bins
        totnevt += binnedNevtRnd

        if verbose:
            print "timebin       = %s-%s ms" % (boundsMin, boundsMax)
            print "Nevt (theor.) =", binnedNevt
            print "Nevt (actual) =", binnedNevtRnd
            print "**************************************"

        if binnedNevtRnd == 0: continue

        # create binned values for energy, mean squared energy and shape parameter
        binnedEnergy = interpolatedEnergy(boundsMin + 0.5*binWidth)
        binnedMSEnergy = interpolatedMSEnergy(boundsMin + 0.5*binWidth)
        binnedAlpha = (2*binnedEnergy**2 - binnedMSEnergy)/(binnedMSEnergy - binnedEnergy**2)

        # define particle for each event in time interval
        for i in range(binnedNevtRnd):
            # Define properties of the particle
            t = boundsMin + random.random() * binWidth
            eNu = get_eNu(binnedAlpha, binnedEnergy)
            (dirx, diry, dirz) = get_direction(eNu)
            ene = get_eE(eNu, dirz)
            # write [t, pid, energy, dirx, diry, dirz] out to file
            outfile.write("%f, %d, %f, %f, %f, %f\n" % (t, pid, ene, dirx, diry, dirz))

    print(("Wrote %i particles to " % totnevt) + output)

    outfile.close()
