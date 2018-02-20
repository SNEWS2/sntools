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
    dSigma_dE = getattr(channel_module, "dSigma_dE")
    get_eE = getattr(channel_module, "get_eE")
    dSigma_dCosT = getattr(channel_module, "dSigma_dCosT")
    bounds_eE = getattr(channel_module, "bounds_eE")
    bounds_eNu = getattr(channel_module, "bounds_eNu")
    targets_per_molecule = getattr(channel_module, "targets_per_molecule")
    pid = getattr(channel_module, "pid")


    '''
    Astrophysics section.
    * neutrinos are well described by a Gamma distribution (arXiv:1211.3920)
    * calculate energy-dependent flux at a fiducial distance of 10 kpc
    '''
    # Convert fiducial distance of 10 kpc into units of MeV**(-1)
    # see http://www.wolframalpha.com/input/?i=10+kpc%2F(hbar+*+c)+in+MeV%5E(-1)
    dSquared = (1.563738e+33)**2

    def dFlux_dE(eNu, time):
        (a, e_sq, luminosity) = flux[time]
        alpha = (2*a**2 - e_sq) / (e_sq - a**2)
        # energy of neutrinos follows a gamma distribution
        gamma_dist = eNu**alpha / gamma(alpha + 1) * ((alpha + 1)/a)**(alpha + 1) * exp(-(alpha + 1)* eNu/a)

        # The `normalization` factor takes into account the oscillation probability
        # as well as the distance (if not equal to 10 kpc).
        return 1/(4*pi*dSquared) * luminosity/a * gamma_dist * normalization


    '''
    Preparation section.
    * Parse input data.
    * For each time step in the input data, calculate instantaneous event rate.
    * Interpolate to get event rate as a function of time.
    '''
    # double differential event rate
    def ddEventRate(eE, eNu, time):
        return dSigma_dE(eNu, eE) * dFlux_dE(eNu, time)


    flux = {}
    nevtValues = []
    molecules_per_kt = 3.343e+31 # number of water molecules in one kt (assuming 18 g/mol)
    n_targets = targets_per_molecule * molecules_per_kt * detectors[detector]

    for (t, mean_e, mean_e_sq, lum) in indata:
        # get time, mean energy, mean squared energy, luminosity
        t = 1000 * t # convert to ms
        flux[t] = (mean_e, mean_e_sq, lum * 624.151) # convert lum from erg/s to MeV/ms

        # integrate over eE and then eNu to obtain the event rate at time t
        simnevt = n_targets * integrate.nquad(ddEventRate, [bounds_eE, bounds_eNu], args=[t]) [0]

        # create a list of nevt values at time t for input into interpolation function
        nevtValues.append(simnevt)

    _flux = sorted([(k,)+v for (k,v) in flux.items()]) # list of tuples: (t, e, e_sq, lum)
    (raw_t, raw_e, raw_e_sq) = [[entry[i] for entry in _flux] for i in range(3)]

    # interpolate the mean energy, mean squared energy and event rate
    interpolatedEnergy = interpolate.pchip(raw_t, raw_e)
    interpolatedMSEnergy = interpolate.pchip(raw_t, raw_e_sq)
    interpolatedNevt = interpolate.pchip(raw_t, nevtValues)


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
    def get_eNu(time):
        dist = lambda _eNu: integrate.quad(ddEventRate, *bounds_eE(_eNu), args=(_eNu, time))[0]
        eNu = rejection_sample(dist, *bounds_eNu, n_bins=200)
        return eNu

    # return direction of scattered electron, if incoming neutrino moves in z direction
    def get_direction(eNu):
        dist = lambda x: dSigma_dCosT(eNu, x)
        cosT = rejection_sample(dist, -1, 1, 200)
        sinT = sin(acos(cosT))
        phi = 2 * pi * random.random() # randomly distributed in [0, 2 pi)
        return (sinT*cos(phi), sinT*sin(phi), cosT)

    bin_width = 1 # bin width in ms
    n_bins = int(duration/bin_width) # number of full-width bins; int() implies floor()
    if verbose:
        print "Now generating events in", bin_width, "ms bins from", starttime, "to", endtime, "ms"
        print "**************************************"

    # scipy is optimized for parallel operation on large arrays, making
    # it orders of magnitude faster to evaluate these interpolated
    # functions for all bins at the same time.
    binned_t = [starttime + (i+0.5)*bin_width for i in range(n_bins)]
    binned_nevt_th = interpolatedNevt(binned_t)
    binned_nevt = np.random.poisson(binned_nevt_th) # Get random number of events in each bin from Poisson distribution
    binned_e = interpolatedEnergy(binned_t)
    binned_e_sq = interpolatedMSEnergy(binned_t)

    for (t, mean_e, mean_e_sq) in zip(binned_t, binned_e, binned_e_sq):
        # we don't need to calculate the luminosity -- it's just a constant
        # factor that doesn't change the result of rejection sampling below
        flux[t] = (mean_e, mean_e_sq, 1)

    with open(output, 'w') as outfile:
        # Iterate over bins to generate events.
        for i in range(n_bins):
            boundsMin = starttime + i * bin_width
            boundsMax = starttime + (i+1) * bin_width

            if verbose:
                print "timebin       = %s-%s ms" % (boundsMin, boundsMax)
                print "Nevt (theor.) =", binned_nevt_th[i]
                print "Nevt (actual) =", binned_nevt[i]
                print "**************************************"

            # define particle for each event in time interval
            for _ in range(binned_nevt[i]):
                # Define properties of the particle
                t = boundsMin + random.random() * bin_width
                eNu = get_eNu(binned_t[i])
                (dirx, diry, dirz) = get_direction(eNu)
                ene = get_eE(eNu, dirz)
                # write [t, pid, energy, dirx, diry, dirz] out to file
                outfile.write("%f, %d, %f, %f, %f, %f\n" % (t, pid, ene, dirx, diry, dirz))

    print "Wrote", sum(binned_nevt), "particles to", output
