"""Parse Garching fluxes.

Input file contains data in the format
  time (s), mean energy (MeV), mean squared energy (MeV^2), luminosity (erg/s)
with one entry per line. We assume that the data follows a gamma distribution
(see arXiv:1211.3920), which is fully described by these values.
See the file 'sample-in.txt' for details.
"""

from math import ceil, floor, gamma, exp
from scipy import interpolate


def parse_input(input, starttime=None, endtime=None):
    """Read simulations data from input file.

    Arguments:
    input -- name of file containing neutrino fluxes
    starttime -- start time set by user via command line option (or None)
    endtime -- end time set by user via command line option (or None)
    """

    # Ensure changes to these variables endure beyond this function's scope.
    global flux, interpolated_e, interpolated_e_sq

    # read data from input file, ignoring lines with comments and empty lines
    with open(input) as infile:
        raw_indata = [map(float, line.split(",")) for line in infile if not (line.startswith("#") or line.isspace())]

    # Compare start/end time entered by user with first/last line of input file
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

    # Ignore data outside of the requested time span.
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

    flux = {}
    for (t, mean_e, mean_e_sq, lum) in indata:
        # save mean energy, mean squared energy, luminosity to dictionary to look up in nu_emission() below
        t = 1000 * t # convert to ms
        flux[t] = (mean_e, mean_e_sq, lum * 624.151) # convert lum from erg/s to MeV/ms

    _flux = sorted([(k,)+v for (k,v) in flux.items()]) # list of tuples: (t, e, e_sq, lum)
    (raw_t, raw_e, raw_e_sq) = [[entry[i] for entry in _flux] for i in range(3)]

    # interpolate the mean energy and mean squared energy
    interpolated_e = interpolate.pchip(raw_t, raw_e)
    interpolated_e_sq = interpolate.pchip(raw_t, raw_e_sq)

    return (starttime, endtime, raw_t)


def prepare_evt_gen(binned_t):
    """Pre-compute values necessary for event generation.

    Scipy/numpy are optimized for parallel operation on large arrays, making
    it orders of magnitude faster to pre-compute all values at one time
    instead of computing them lazily when needed.

    Argument:
    binned_t -- list of time bins for generating events
    """
    binned_e = interpolated_e(binned_t)
    binned_e_sq = interpolated_e_sq(binned_t)

    for (t, mean_e, mean_e_sq) in zip(binned_t, binned_e, binned_e_sq):
        # Since the event rate is already known, luminosity can now be set to a
        # constant. We save time by not calculating it.
        flux[t] = (mean_e, mean_e_sq, 1)

    return None


def nu_emission(eNu, time):
    """Number of neutrinos emitted, as a function of energy.

    This is not yet the flux! The geometry factor 1/(4 pi r**2) is added later.
    Arguments:
    eNu -- neutrino energy
    time -- time ;)
    """
    (e, e_sq, luminosity) = flux[time]
    alpha = (2 * e**2 - e_sq) / (e_sq - e**2)

    # total number = luminosity / mean energy
    number = luminosity / e
    # energy of neutrinos follows a gamma distribution
    gamma_dist = eNu**alpha / gamma(alpha + 1) * ((alpha + 1)/e)**(alpha + 1) * exp(-(alpha + 1) * eNu/e)

    return number * gamma_dist
