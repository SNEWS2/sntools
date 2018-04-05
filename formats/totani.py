"""Parse Totani fluxes.

Based on simulation from arXiv:astro-ph/9710203, which is the primary
model used in the Hyper-Kamiokande Design Report.
"""

from math import ceil, floor, gamma, exp, log10
from scipy import interpolate

zero = 1E-99 # not exactly zero to ensure log interpolation is still possible


def parse_input(input, inflv, starttime=None, endtime=None):
    """Read simulations data from input file.

    Arguments:
    input -- prefix of file containing neutrino fluxes
    inflv -- neutrino flavor to consider
    starttime -- start time set by user via command line option (or None)
    endtime -- end time set by user via command line option (or None)
    """
    global times, e_bins
    global N_dict, egroup_dict, dNLde_dict, log_spectrum
    times = []
    e_bins = [zero] # energy bins are the same for all times; first bin = 0 MeV
    N_dict, egroup_dict, dNLde_dict, log_spectrum = {}, {}, {}, {}

    # _parse() is a helper function to parse the files since their format is
    # not straightforward. It fills in the global variables defined above.
#     if inflv == "e":
#         # nu_e fluxes during the neutronization burst are in a separate file,
#         # with a somewhat different format
#         _parse_nb(input + "-nb.txt")
#
    _parse(input + "-early.txt", "early", inflv)
    _parse(input + "-late.txt", "late", inflv)

    # Compare start/end time entered by user with first/last line of input file
    _starttime = times[0]
    _endtime = times[-1]

    if not starttime:
        starttime = ceil(_starttime)
    elif starttime < _starttime:
        print("Error: Start time must be greater than earliest time in input files. Aborting ...")
        exit()

    if not endtime:
        endtime = floor(_endtime)
    elif endtime > _endtime:
        print("Error: End time must be less than latest time in input files. Aborting ...")
        exit()

    # if user entered a custom start/end time, find indices of relevant time bins
    i_min, i_max = 0, len(times) - 1
    for (i, time) in enumerate(times):
        if time < starttime:
            i_min = i
        elif time > endtime:
            i_max = i
            break

    # fill in global dictionaries that will get used later
    _additional_setup()

    return (starttime, endtime, times[i_min:i_max+1])


def prepare_evt_gen(binned_t):
    """Pre-compute values necessary for event generation.

    Scipy/numpy are optimized for parallel operation on large arrays, making
    it orders of magnitude faster to pre-compute all values at one time
    instead of computing them lazily when needed.

    Argument:
    binned_t -- list of time bins for generating events
    """
    for time in binned_t:
        if dNLde_dict.has_key(time):
            # we have already computed dNLde and the interpolated spectrum
            # for this time bin in _additional_setup() below
            continue

        # find closest time bins -> t0, t1
        for t_bin in times:
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
            tmp = prev_dNLde[i] + (next_dNLde[i] - prev_dNLde[i]) * (time-t0)/(t1-t0)
            dNLde.append(tmp)

        dNLde_dict[time] = dNLde

        # Get emission spectrum by log cubic spline interpolation
        log_group_e = [log10(e_bin) for e_bin in e_bins]
        log_dNLde = [log10(d) for d in dNLde]
        log_spectrum[time] = interpolate.pchip(log_group_e, log_dNLde)

    return None

def nu_emission(eNu, time):
    """Number of neutrinos emitted, as a function of energy.

    This is not yet the flux! The geometry factor 1/(4 pi r**2) is added later.
    Arguments:
    eNu -- neutrino energy
    time -- time ;)
    """
    f = log_spectrum[time]
    return 10 ** f(log10(eNu)) # transform log back to actual value


"""Helper functions.

(Which I'll need. Because the file format is ... not exactly pretty. -.-)
"""
def _parse(input, format, flv):
    """Read data from files into dictionaries to look up by time."""
    with open(input) as infile:
        raw_indata = [line for line in infile]

    chunks = []

    if format == "early":
        # 42 lines per time bin, 26 bins in wilson-early.txt
        for i in range(26):
            chunks.append(raw_indata[42*i:42*(i+1)])
        line_N = 6
        range_egroup = range(19, 39)
    elif format == "late":
        # 46 lines per time bin, 36 bins in wilson-late.txt
        for i in range(36):
            chunks.append(raw_indata[46*i:46*(i+1)])
        line_N = 8
        range_egroup = range(21, 41)

    # input files contain information for e, eb & x right next to each other,
    # so depending on the flavor, we might need an offset
    offset = {"e": 0, "eb": 1, "x": 2, "xb": 2}[flv]

    # for each time bin, save data to dictionaries to look up later
    for chunk in chunks:
        # first line contains time
        time = float(chunk[0].split()[0]) * 1000 # convert to ms
        times.append(time)

        # N = total number of neutrinos emitted up to this time
        N = float(chunk[line_N].split()[offset])
        if offset == 2: N /= 4 # file contains sum of nu_mu, nu_tau and anti-particles
        N_dict[time] = N

        # number of neutrinos emitted in this time bin, separated into energy bins
        egroup = [zero] # start with 0 neutrinos at 0 MeV bin
        for i in range_egroup:
            line = map(float, chunk[i].split())
            egroup.append(line[-3+offset])

            # Once, for the very first time bin, save the energy bins:
            if egroup_dict == {}:
                e_bins.append(line[1] / 1000) # energy of this bin (in MeV)

        egroup_dict[time] = egroup

    return None


def _parse_nb(input):
    # TODO: read in data from wilson-nb.txt and fill dictionaries
    # Note: that file comes from a different simulation, so ...
    #       * the simulation time is inconsistent -> use time after bounce!
    #       * there's a discontinuity between it and data in the early/late file,
    #         so I have to interpolate (see README.txt)
    return None


def _additional_setup():
    """Calculate energy spectrum and number luminosity for each time bin."""
    for (i, time) in enumerate(times):
        # Get energy spectrum per MeV^-1 instead of in (varying-size) energy bins
        E_integ = 0
        spec = []
        egroup = egroup_dict[time] # list: number of neutrinos in different e_bins

        for (j, n) in enumerate(egroup):
            if j == 0 or j == len(egroup)-1:
                spec.append(zero)
            else:
                spec.append(n / (e_bins[j+1] - e_bins[j-1]))
                E_integ += (spec[j-1] + spec[j]) * (e_bins[j] - e_bins[j-1]) / 2

        spec = [x / E_integ for x in spec] # normalise to 1

        # Calculate number luminosity
        if i == 0:
            # TODO: How does this interact with data from wilson-nb.txt?
            num_lum = zero
        else:
            prev_time = times[i-1]
            num_lum = (N_dict[time] - N_dict[prev_time]) / (time - prev_time)

        # Calculate differential number luminosity
        dNLde = [num_lum * spectrum for spectrum in spec]
        dNLde_dict[time] = dNLde

        # Get emission spectrum by log cubic spline interpolation
        log_group_e = [log10(e_bin) for e_bin in e_bins]
        log_dNLde = [log10(d) for d in dNLde]
        log_spectrum[time] = interpolate.pchip(log_group_e, log_dNLde)

    return None
