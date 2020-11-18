"""Parse Totani fluxes.

Based on simulation from arXiv:astro-ph/9710203, which is the primary
model used in the Hyper-Kamiokande Design Report.
All times are post-bounce, which differs slightly from the simulation time given
in the raw input files. We fix this when reading from the files in _parse().
"""

from math import log10
from scipy.interpolate import InterpolatedUnivariateSpline
from sntools.formats import get_endtime, get_starttime

zero = 1e-99  # not exactly zero to ensure log interpolation is still possible


def parse_input(input, inflv, starttime, endtime):
    """Read simulations data from input file.

    Arguments:
    input -- prefix of file containing neutrino fluxes
    inflv -- neutrino flavor to consider
    starttime -- start time set by user via command line option (or None)
    endtime -- end time set by user via command line option (or None)
    """
    global times, times_el, times_nb, e_bins
    global N_dict, egroup_dict, dNLde_dict, log_spectrum
    times_el, times_nb = [], []
    e_bins = [zero]  # energy bins are the same for all times; first bin = 0 MeV
    N_dict, egroup_dict, dNLde_dict, log_spectrum = {}, {}, {}, {}

    # The file format is complicated, so we define helper functions below
    _parse(input + "-early.txt", "early", inflv)
    _parse(input + "-late.txt", "late", inflv)
    _calculate_dNLde()  # calculate number luminosity for early and late files
    times = times_el

    starttime = get_starttime(starttime, times[0])
    endtime = get_endtime(endtime, times[-1])

    # If user entered a custom start/end time, select only relevant time bins
    i_min, i_max = 0, len(times) - 1
    for (i, time) in enumerate(times):
        if time < starttime:
            i_min = i
        elif time > endtime:
            i_max = i
            break
    times = times[i_min : i_max + 1]

    # nu_e fluxes during the neutronization burst (40-50 ms) are in a separate
    # file, with more precise time bins and a different format.
    if inflv == "e" and starttime < 50:
        _parse_nb(input + "-nb.txt")
        times = sorted(times + times_nb)

    # Get spectra for relevant time bins by log cubic spline interpolation
    log_group_e = [log10(e_bin) for e_bin in e_bins]
    for time in times:
        log_dNLde = [log10(d) for d in dNLde_dict[time]]
        log_spectrum[time] = InterpolatedUnivariateSpline(log_group_e, log_dNLde)

    return (starttime, endtime, times)


def prepare_evt_gen(binned_t):
    """Pre-compute values necessary for event generation.

    Scipy/numpy are optimized for parallel operation on large arrays, making
    it orders of magnitude faster to pre-compute all values at one time
    instead of computing them lazily when needed.

    Argument:
    binned_t -- list of time bins for generating events
    """
    for time in binned_t:
        if time in log_spectrum:
            # we have already computed the interpolated spectrum at this time
            continue

        if 40 <= time <= 49.99 and times_nb != []:
            # take fluxes from nb file into account
            _times = [x for x in times_nb if x in times]
        else:  # use fluxes from early/late file
            _times = [x for x in times_el if x in times]

        # find closest time bins -> t0, t1
        for t_bin in _times:
            if time <= t_bin:
                t1 = t_bin
                break
            else:
                t0 = t_bin

        # get dNLde at the intermediate time
        dNLde = []
        prev_dNLde = dNLde_dict[t0]
        next_dNLde = dNLde_dict[t1]
        for (i, _) in enumerate(e_bins):
            # linear interpolation over time each energy bin
            tmp = prev_dNLde[i] + (next_dNLde[i] - prev_dNLde[i]) * (time - t0) / (t1 - t0)
            dNLde.append(tmp)

        # Get emission spectrum by log cubic spline interpolation
        log_group_e = [log10(e_bin) for e_bin in e_bins]
        log_dNLde = [log10(d) for d in dNLde]
        log_spectrum[time] = InterpolatedUnivariateSpline(log_group_e, log_dNLde)

    return None


def nu_emission(eNu, time):
    """Number of neutrinos emitted, as a function of energy.

    This is not yet the flux! The geometry factor 1/(4 pi r**2) is added later.
    Arguments:
    eNu -- neutrino energy
    time -- time ;)
    """
    f = log_spectrum[time]
    return 10 ** f(log10(eNu))  # transform log back to actual value


# Helper functions
def _parse(input, format, flv):
    """Read data from files into dictionaries to look up by time."""
    with open(input) as infile:
        raw_indata = [line for line in infile]

    chunks = []

    if format == "early":
        # 42 lines per time bin, 26 bins in wilson-early.txt
        for i in range(26):
            chunks.append(raw_indata[42 * i : 42 * (i + 1)])
        line_N = 6
        range_egroup = range(19, 39)
    elif format == "late":
        # 46 lines per time bin, 36 bins in wilson-late.txt
        for i in range(36):
            chunks.append(raw_indata[46 * i : 46 * (i + 1)])
        line_N = 8
        range_egroup = range(21, 41)

    # input files contain information for e, eb & x right next to each other,
    # so depending on the flavor, we might need an offset
    offset = {"e": 0, "eb": 1, "x": 2, "xb": 2}[flv]

    # for each time bin, save data to dictionaries to look up later
    for chunk in chunks:
        # first line contains time
        time = float(chunk[0].split()[0]) * 1000  # convert to ms
        time -= 2  # change from simulation time into time after core bounce
        times_el.append(time)

        # N = total number of neutrinos emitted up to this time
        N = float(chunk[line_N].split()[offset])
        if offset == 2:
            N /= 4  # file contains sum of nu_mu, nu_tau and anti-particles
        N_dict[time] = N

        # number of neutrinos emitted in this time bin, separated into energy bins
        egroup = [zero]  # start with 0 neutrinos at 0 MeV bin
        for i in range_egroup:
            line = list(map(float, chunk[i].split()))
            egroup.append(line[-3 + offset])

            # Once, for the very first time bin, save the energy bins:
            if egroup_dict == {}:
                e_bins.append(line[1] / 1000)  # energy of this bin (in MeV)

        egroup_dict[time] = egroup

    return None


def _parse_nb(input):
    """More granular nu_e data for the neutronization burst ("nb", 40-50ms).

    Note: the nb file comes from a slightly different simulation, therefore we
    have to deal with a time offset and a scaling factor.
    """
    with open(input) as infile:
        raw_indata = [line for line in infile]

    # 26 lines per time bin, 99 bins in wilson-nb.txt. Bin 6 is equivalent to
    # 40ms post-bounce & bin 56 is 50ms, so we only select that range:
    chunks = [raw_indata[26 * i : 26 * (i + 1)] for i in range(6, 57)]

    # for each time bin, save data to dictionaries to look up later
    for chunk in chunks:
        time = float(chunk[0].split()[2]) * 1000  # convert to ms
        time -= 467.5  # 40-50ms post-bounce equals 507.5-517.5ms in this file
        times_nb.append(time)

        luminosity = float(chunk[1].split()[2]) * 624.151  # convert erg/s to MeV/ms

        # number of neutrinos emitted in this time bin, separated into energy bins
        egroup = [zero]  # start with 0 neutrinos at 0 MeV bin
        for i in range(3, 23):
            line = list(map(float, chunk[i].split()))
            egroup.append(line[-3])

        # Get energy spectrum per MeV^-1 instead of in (varying-size) energy bins
        E_integ = 0
        spec = []
        for (j, n) in enumerate(egroup):
            if j == 0 or j == len(egroup) - 1:
                spec.append(zero)
            else:
                spec.append(n / (e_bins[j + 1] - e_bins[j - 1]))
                E_integ += (spec[j - 1] * e_bins[j - 1] + spec[j] * e_bins[j]) * (e_bins[j] - e_bins[j - 1]) / 2

        spec = [x / E_integ * luminosity for x in spec]

        # nb and early/late data come from slightly different simulations and
        # have a discontinuity, so we scale with a time-dependent factor
        nb_scale = 1 - 5.23 / 13.82 * (time - 40) / 10  # 1 at 40ms, 8.59/13.82 at 50ms
        dNLde_dict[time] = [x * nb_scale for x in spec]
    return None


def _calculate_dNLde():
    """Calculate number luminosity spectrum for each time bin."""
    for (i, time) in enumerate(times_el):
        # Get energy spectrum per MeV^-1 instead of in (varying-size) energy bins
        E_integ = 0
        spec = []
        egroup = egroup_dict[time]  # list: number of neutrinos in different e_bins

        for (j, n) in enumerate(egroup):
            if j == 0 or j == len(egroup) - 1:
                spec.append(zero)
            else:
                spec.append(n / (e_bins[j + 1] - e_bins[j - 1]))
                E_integ += (spec[j - 1] + spec[j]) * (e_bins[j] - e_bins[j - 1]) / 2

        spec = [x / E_integ for x in spec]  # normalise to 1

        # Calculate number luminosity
        if i == 0:
            num_lum = zero
        else:
            prev_time = times_el[i - 1]
            num_lum = (N_dict[time] - N_dict[prev_time]) / (time - prev_time)

        # Calculate differential number luminosity
        dNLde = [num_lum * spectrum for spectrum in spec]
        dNLde_dict[time] = dNLde

    return None
