#!/usr/bin/python

from importlib import import_module
from math import pi, sin, cos, acos
import numpy as np
import random
from scipy import integrate, interpolate

def gen_evts(channel, input, format, inflv, normalization, starttime, endtime, verbose):

    """Setup.

    * Import module for reading input file.
    * Import module with interaction channel-specific functions & parameters.
    * Other setup.
    """

    format = import_module("formats." + format)
    channel = import_module("interaction_channels." + channel)

    # The `normalization` factor includes oscillation probability, distance of
    # supernova (if not equal to 10 kpc) and number of water molecules in detector.
    normalization *= channel.targets_per_molecule


    """Helper functions."""
    # double differential event rate
    def ddEventRate(eE, eNu, time):
        return channel.dSigma_dE(eNu, eE) * dFlux_dE(eNu, time)

    # dFlux_dE(eNu, time) is called hundreds of times for each generated event,
    # often with repetitive arguments (when integrating ddEventRate over eE).
    # To save time, we cache results in a dictionary.
    cached_flux = {}
    def dFlux_dE(eNu, time):
        if not cached_flux.has_key((eNu, time)):
            fiducial_distance = 1.563738e+33 # 10 kpc/(hbar * c) in MeV**(-1)
            emission = format.nu_emission(eNu, time)
            cached_flux[(eNu, time)] = emission / (4 * pi * fiducial_distance**2)
        return cached_flux[(eNu, time)]

    # get a value from an arbitrary distribution dist
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

    # use rejection sampling to get the energy of an interacting neutrino
    def get_eNu(time):
        dist = lambda _eNu: integrate.quad(ddEventRate, *channel.bounds_eE(_eNu), args=(_eNu, time))[0]
        eNu = rejection_sample(dist, *channel.bounds_eNu, n_bins=200)
        return eNu

    # get direction of outgoing particle (incoming neutrino moves in z direction)
    def get_direction(eNu):
        dist = lambda _cosT: channel.dSigma_dCosT(eNu, _cosT)
        cosT = rejection_sample(dist, -1, 1, 200)
        sinT = sin(acos(cosT))
        phi = 2 * pi * random.random() # randomly distributed in [0, 2 pi)
        return (sinT*cos(phi), sinT*sin(phi), cosT)


    """Main section.

    * Get event rate by interpolating from time steps in the input data.
    * For each 1ms bin, get number of events from a Poisson distribution.
    * Generate these events from time-dependent energy & direction distribution.
    """
    (starttime, endtime, raw_times) = format.parse_input(input, inflv, starttime, endtime)

    # integrate over eE and then eNu to obtain the event rate at time t
    raw_nevts = [normalization * integrate.nquad(ddEventRate, [channel.bounds_eE, channel.bounds_eNu], args=[t])[0]
                 for t in raw_times]
    event_rate = interpolate.pchip(raw_times, raw_nevts)

    bin_width = 1 # in ms
    n_bins = int((endtime - starttime)/bin_width) # number of full-width bins; int() implies floor()
    if verbose:
        print "Now generating events in", bin_width, "ms bins from", starttime, "to", endtime, "ms"
        print "**************************************"

    # scipy is optimized for operating on large arrays, making it orders of
    # magnitude faster to pre-compute all values of the interpolated functions.
    binned_t = [starttime + (i+0.5)*bin_width for i in range(n_bins)]
    binned_nevt_th = event_rate(binned_t)
    binned_nevt = np.random.poisson(binned_nevt_th) # Get random number of events in each bin from Poisson distribution
    format.prepare_evt_gen(binned_t) # give flux script a chance to pre-compute values

    events = []
    for i in range(n_bins):
        t0 = starttime + i * bin_width

        if verbose and i%(10**(4-verbose)) == 0:
            print "%s-%s ms: %d events (%.5f expected)" % (t0, t0+bin_width, binned_nevt[i], binned_nevt_th[i])

        # generate events in this time bin
        for _ in range(binned_nevt[i]):
            t = t0 + random.random() * bin_width
            eNu = get_eNu(binned_t[i])
            (dirx, diry, dirz) = get_direction(eNu)
            ene = channel.get_eE(eNu, dirz)
            if verbose and ene < thr_e: thr_nevt -= 1
            events.append((t, channel.pid, ene, dirx, diry, dirz))

    print "Generated %s particles (expected: %.2f particles)" % (sum(binned_nevt), sum(binned_nevt_th))

    return events
