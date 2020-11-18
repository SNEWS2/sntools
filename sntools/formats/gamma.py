"""Parse Gamma fluxes.

Input file contains data in the format
  time (s), mean energy (MeV), mean squared energy (MeV^2), luminosity (erg/s)
with one entry per line. We assume that the data follows a gamma distribution
(see arXiv:1211.3920), which is fully described by these values.
See the file 'fluxes/sample-gamma.txt' for details.
"""

from math import gamma, exp
from scipy import interpolate
from sntools.formats import get_endtime, get_starttime


def parse_input(input, inflv, starttime, endtime):
    """Read simulations data from input file.

    Arguments:
    input -- prefix of file containing neutrino fluxes
    inflv -- neutrino flavor to consider
    starttime -- start time set by user via command line option (or None)
    endtime -- end time set by user via command line option (or None)
    """
    # read data from input file, ignoring lines with comments and empty lines
    with open(input) as infile:
        raw_indata = [
            list(map(float, line.split(","))) for line in infile if not (line.startswith("#") or line.isspace())
        ]
    for entry in raw_indata:
        entry[0] *= 1000  # convert time to ms

    starttime = get_starttime(starttime, raw_indata[0][0])
    endtime = get_endtime(endtime, raw_indata[-1][0])

    # Ignore data outside of the requested time span.
    indata = []
    for (i, entry) in enumerate(raw_indata):
        if i == 0:
            continue
        if entry[0] > starttime:
            indata.append(raw_indata[i - 1])
            if entry[0] > endtime:
                indata.append(entry)
                break

    # save mean energy, mean squared energy, luminosity to dictionary to look up in nu_emission() below
    global flux
    flux = {}
    for timebin in indata:
        # input files contain information for nu_e in columns 1-3, for
        # anti-nu_e in cols 4-6 and for nu_x in columns 7-9
        offset = {"e": 1, "eb": 4, "x": 7, "xb": 7}[inflv]
        (mean_e, mean_e_sq, lum) = timebin[offset : offset + 3]
        t = timebin[0]
        flux[t] = (mean_e, mean_e_sq, lum * 624.151)  # convert lum from erg/s to MeV/ms

    return (starttime, endtime, sorted(flux.keys()))


def prepare_evt_gen(binned_t):
    """Pre-compute values necessary for event generation.

    Scipy/numpy are optimized for parallel operation on large arrays, making
    it orders of magnitude faster to pre-compute all values at one time
    instead of computing them lazily when needed.

    Argument:
    binned_t -- list of time bins for generating events
    """
    _flux = sorted([(k,) + v for (k, v) in flux.items()])  # list of tuples: (t, e, e_sq, lum)
    (raw_t, raw_e, raw_e_sq, raw_lum) = [[entry[i] for entry in _flux] for i in range(4)]

    # interpolate mean energy, mean squared energy and luminosity ...
    interpolated_e = interpolate.pchip(raw_t, raw_e)
    interpolated_e_sq = interpolate.pchip(raw_t, raw_e_sq)
    interpolated_lum = interpolate.pchip(raw_t, raw_lum)
    # ... and evaluate them at all relevant times
    binned_e = interpolated_e(binned_t)
    binned_e_sq = interpolated_e_sq(binned_t)
    binned_lum = interpolated_lum(binned_t)

    for (t, mean_e, mean_e_sq, mean_lum) in zip(binned_t, binned_e, binned_e_sq, binned_lum):
        flux[t] = (mean_e, mean_e_sq, mean_lum)

    return None


def nu_emission(eNu, time):
    """Number of neutrinos emitted, as a function of energy.

    This is not yet the flux! The geometry factor 1/(4 pi r**2) is added later.
    Arguments:
    eNu -- neutrino energy
    time -- time ;)
    """
    (e, e_sq, luminosity) = flux[time]
    alpha = (2 * e ** 2 - e_sq) / (e_sq - e ** 2)

    # energy of neutrinos follows a gamma distribution
    gamma_dist = eNu ** alpha / gamma(alpha + 1) * ((alpha + 1) / e) ** (alpha + 1) * exp(-(alpha + 1) * eNu / e)
    # total number = luminosity / mean energy
    return luminosity / e * gamma_dist
