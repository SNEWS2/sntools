#!/usr/bin/python

from math import pi, sin, cos, acos
import numpy as np
import random
from scipy import integrate, interpolate


def gen_evts(_channel, _flux, n_targets, seed, verbose):
    """Generate events.

    * Get event rate by interpolating from time steps in the input data.
    * For each 1ms bin, get number of events from a Poisson distribution.
    * Generate these events from time-dependent energy & direction distribution.

    Arguments:
    _channel -- BaseChannel instance for the current interaction channel
    _flux -- BaseFlux instance with appropriate flavor and time range (includes weighting due to flux transformation and distance)
    n_targets -- number of target particles in detector
    seed -- random number seed to reproducibly generate events
    """
    random.seed(seed)
    np.random.seed(int(seed))

    global channel, cached_flux, flux
    flux = _flux
    channel = _channel
    tag = str(channel.__class__).split('.')[-2]
    if tag in ('c12nc', 'es'):
        tag += '-' + str(channel).split("'")[-2]

    # ddEventRate(eE, eNu, time) is called hundreds of times for each generated event,
    # often with identical eNu and time values (when integrating over eE).
    # To save time, we cache results in a dictionary.
    cached_flux = {}

    # integrate over eE and then eNu to obtain the event rate at time t
    if verbose:
        print(f"[{tag}] Calculating event rate for {flux} ...")
    raw_nevts = [n_targets * integrate.nquad(ddEventRate, [channel.bounds_eE, channel.bounds_eNu], args=[t], opts=[channel._opts, {}])[0]
                 for t in flux.raw_times]
    event_rate = interpolate.pchip(flux.raw_times, raw_nevts)

    bin_width = 1  # in ms
    n_bins = int((flux.endtime - flux.starttime) / bin_width)  # number of full-width bins; int() implies floor()
    if verbose:
        print(f"[{tag}] Generating events in {bin_width} ms bins from {flux.starttime} to {flux.endtime} ms ...")

    # scipy is optimized for operating on large arrays, making it orders of
    # magnitude faster to pre-compute all values of the interpolated functions.
    binned_t = [flux.starttime + (i + 0.5) * bin_width for i in range(n_bins)]
    binned_nevt_th = event_rate(binned_t)
    # check for unphysical values of interpolated function event_rate(t)
    for _i, _n in enumerate(binned_nevt_th):
        if _n < 0:
            binned_nevt_th[_i] = 0
    binned_nevt = np.random.poisson(binned_nevt_th)  # Get random number of events in each bin from Poisson distribution
    flux.prepare_evt_gen(binned_t)  # give flux script a chance to pre-compute values

    events = []
    for i in range(n_bins):
        t0 = flux.starttime + i * bin_width

        if verbose and i % (10 ** (4 - verbose)) == 0:
            print(f"[{tag}] {t0}-{t0 + bin_width} ms: {binned_nevt[i]} events ({binned_nevt_th[i]:.5f} expected)")

        # generate events in this time bin
        for _ in range(binned_nevt[i]):
            eNu = get_eNu(binned_t[i])
            direction = get_direction(eNu)  # (dirx, diry, dirz)
            evt = channel.generate_event(eNu, *direction)
            evt.time = t0 + random.random() * bin_width
            events.append(evt)

    print(f"[{tag}] Generated {sum(binned_nevt)} particles (expected: {sum(binned_nevt_th):.2f} particles)")

    return events


# Helper functions
def ddEventRate(eE, eNu, time):
    """Double differential event rate."""
    if (eNu, time) not in cached_flux:
        cached_flux[(eNu, time)] = flux.nu_emission(eNu, time)
    return channel.dSigma_dE(eNu, eE) * cached_flux[(eNu, time)]


def rejection_sample(dist, min_val, max_val, n_bins=100):
    """Sample value from an arbitrary distribution."""
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
    for j in range(max(j_max - 9, 0), min(j_max + 10, n_bins)):
        val = min_val + bin_width * (j + 0.5)
        p = dist(val)
        if p > p_max:
            p_max = p

    while True:
        val = min_val + (max_val - min_val) * random.random()
        if p_max * random.random() < dist(val):
            break

    return val


def get_eNu(time):
    """Get energy of interacting neutrino using rejection sampling."""
    def dist(eNu):
        return integrate.quad(ddEventRate, *channel.bounds_eE(eNu), args=(eNu, time), points=channel._opts(eNu)["points"])[0]
    eNu = rejection_sample(dist, *channel.bounds_eNu, n_bins=200)
    return eNu


def get_direction(eNu):
    """Get direction of outgoing particle using rejection sampling.
    (Assumes that incoming neutrino with energy eNu moves in z direction.)
    """
    def dist(cosT):
        return channel.dSigma_dCosT(eNu, cosT)
    cosT = rejection_sample(dist, -1, 1, 200)
    sinT = sin(acos(cosT))
    phi = 2 * pi * random.random()  # randomly distributed in [0, 2 pi)
    return (sinT * cos(phi), sinT * sin(phi), cosT)
