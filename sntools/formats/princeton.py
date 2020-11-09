"""Parse Princeton fluxes.

For simulations by Seadrow et al., arXiv:1804.00689
and later re-runs by David Vartanyan.
Format:
time, dL/dE (nu_e), dL/dE (anti-nu_e), dL/dE (nu_x)
where the spectral luminosity is split in 20 energy bins per flavor.
See parsing code below for details.
"""

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
    global times, dNLdE
    times = []
    dNLdE = {}

    with open(input) as infile:
        indata = [list(map(float, line.split())) for line in infile if not (line.startswith("#") or line.isspace())]

    # input files contain information for e, eb & x in neighbouring columns,
    # so depending on the flavor, we might need an offset
    offset = {"e": 1, "eb": 21, "x": 41, "xb": 41}[inflv]

    # luminosity is in 20 bins covering 1-300 MeV (for e), 1-100 MeV (for eb & x)
    emax = 300 if inflv == "e" else 100
    ebins = [0] + [emax ** ((i + 0.5) * 0.05) for i in range(22)]  # add extra bins at start/end for interpolation

    # for each time bin, save data to dictionaries to look up later
    for line in indata:
        time = line[0] * 1000  # convert time to ms
        time -= 31.7  # offset between time in file and core bounce (D. Vartanyan, private communications)
        times.append(time)

        diff_number_flux = [0]  # Set flux at 0 MeV to 0
        for emean, diff_lum in zip(ebins[1:-1], line[offset : offset + 20]):
            diff_lum *= 1e50  # file gives spectral luminosity in 10^50 erg/s/MeV
            diff_lum *= 624.151  # convert erg/s/MeV to MeV/ms/MeV
            if offset == 41:
                diff_lum /= 4  # file contains sum of nu_mu, nu_tau and anti-particles
            number_flux = diff_lum / emean
            diff_number_flux.append(number_flux)
        # Let flux at >100 MeV smoothly go to zero
        diff_number_flux.append(diff_number_flux[-1] * 0.001)
        diff_number_flux.append(0)

        dNLdE[time] = interpolate.pchip(ebins, diff_number_flux)

    starttime = get_starttime(starttime, times[0])
    endtime = get_endtime(endtime, times[-1])

    # if user entered a custom start/end time, find indices of relevant time bins
    i_min, i_max = 0, len(times) - 1
    for (i, time) in enumerate(times):
        if time < starttime:
            i_min = i
        elif time > endtime:
            i_max = i
            break

    return (starttime, endtime, times[i_min : i_max + 1])


def prepare_evt_gen(binned_t):
    """Pre-compute values necessary for event generation.

    Scipy/numpy are optimized for parallel operation on large arrays, making
    it orders of magnitude faster to pre-compute all values at one time
    instead of computing them lazily when needed.

    Argument:
    binned_t -- list of time bins for generating events
    """
    # unnecessary here; linear interpolation is fast enough to do it on demand
    return None


def nu_emission(eNu, time):
    """Number of neutrinos emitted, as a function of energy.

    This is not yet the flux! The geometry factor 1/(4 pi r**2) is added later.
    Arguments:
    eNu -- neutrino energy
    time -- time ;)
    """
    # find previous/next time bin and perform linear interpolation
    for t_prev, t_next in zip(times[:-1], times[1:]):
        if time < t_next:
            break

    dNLdE_prev = dNLdE[t_prev](eNu)
    dNLdE_next = dNLdE[t_next](eNu)
    result = dNLdE_prev + (dNLdE_next - dNLdE_prev) * (time - t_prev) / (t_next - t_prev)

    return result
