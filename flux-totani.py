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
    global E_dict, N_dict, ave_dict, egroup_dict, spec_dict, NL_dict, dNLde_dict, interpolated_spec_dict
    times = []
    e_bins = [zero] # energy bins are the same for all times; first bin = 0 MeV
    # TODO: Some of these dictionaries probably don't need to be global.
    E_dict = {} # this information is never used
    N_dict = {}
    ave_dict = {} # never used
    egroup_dict = {}
    spec_dict = {}
    NL_dict = {}
    dNLde_dict = {}
    interpolated_spec_dict = {}

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
        print("Error: Start time must be greater than earliest time in input file. Aborting ...")
        exit()

    if not endtime:
        endtime = floor(_endtime)
    elif endtime > _endtime:
        print("Error: End time must be less than latest time in input file. Aborting ...")
        exit()

    # fill in global dictionaries that will get used later
    _additional_setup()

    # if user entered a custom start/end time, consider only relevant time bins
    i_min, i_max = 0, len(times)
    for (i, time) in enumerate(times):
        if time < starttime:
            i_min = i
        elif time > endtime:
            i_max = i+1
            break
    times = times[i_min:i_max]

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
        if dNLde_dict.has_key(time):
            # we have already computed dNLde_ and the interpolated spectrum
            # for this time bin in _additional_setup() below
            continue

        # find closest time bins -> t0, t1
        for t_bin in times:
            if time <= t_bin:
                t1 = t_bin
                break
            else:
                t0 = t_bin

        # get dNLde_ at the intermediate time
        dNLde_ = []
        prev_dNLde = dNLde_dict[t0]
        next_dNLde = dNLde_dict[t1]
        for (i, _) in enumerate(e_bins):
            # linear interpolation over time each energy bin
            tmp = prev_dNLde[i] + (next_dNLde[i] - prev_dNLde[i]) * (time-t0)/(t1-t0)
            dNLde_.append(tmp)

        dNLde_dict[time] = dNLde_

        ### Get emission spectrum by log cubic spline interpolation ###
        log_group_e = [log10(e_bin) for e_bin in e_bins] # difference between base 10 and base e is just a constant factor of ln(10)
        log_dNLde_ = [log10(d) for d in dNLde_]
        interpolated_spec_dict[time] = interpolate.pchip(log_group_e, log_dNLde_)

    return None

def nu_emission(eNu, time):
    """Number of neutrinos emitted, as a function of energy.

    This is not yet the flux! The geometry factor 1/(4 pi r**2) is added later.
    Arguments:
    eNu -- neutrino energy
    time -- time ;)
    """
    interpolated_f = interpolated_spec_dict[time]
    return 10 ** interpolated_f(log10(eNu)) # transform log back to actual value


"""Helper functions.

(Which I'll need. Because the file format is ... not exactly pretty. -.-)
"""
def _parse(input, format, flv):
    """Read data from files into dictionaries to look up by time."""
    with open(input) as infile:
        raw_indata = [line for line in infile]

    chunks = []

    if format == "early":
        # 42 lines per chunk of data (i.e. time bin), 26 chunks in wilson-early.txt
        for i in range(26):
            chunks.append(raw_indata[42*i:42*(i+1)])
        line_E = 3
        line_N = 6
        line_ave = 16
        range_egroup = range(19, 39)
    elif format == "late":
        # 46 lines per chunk of data (i.e. time bin), 36 chunks in wilson-late.txt
        for i in range(36):
            chunks.append(raw_indata[46*i:46*(i+1)])
        line_E = 5 # TODO: replace these variables by a constant offset of 2 w/r/t the "early" line numbers?
        line_N = 8
        line_ave = 18
        range_egroup = range(21, 41)

    # input files contain information for e, eb & x right next to each other,
    # so depending on the flavor, we might need an offset
    offset = {"e": 0, "eb": 1, "x": 2, "xb": 2}
    flv = offset[flv]

    # for each time bin, save data to dictionaries to look up later
    for chunk in chunks:
        # first line contains time
        time = float(chunk[0].split()[0]) * 1000 # convert to ms
        times.append(time)

        # E_ = total energy emitted up to this time in erg
        E_ = float(chunk[line_E].split()[flv]) # TODO: convert to MeV?
        if flv == 2: E_ /= 4 # file contains sum of nu_mu, nu_tau and anti-particles
        E_dict[time] = E_

        # N_ = total number of neutrinos emitted up to this time
        N_ = float(chunk[line_N].split()[flv])
        if flv == 2: N_ /= 4 # file contains sum of nu_mu, nu_tau and anti-particles
        N_dict[time] = N_

        # ave_ = average neutrino energy in keV
        ave_ = float(chunk[line_ave].split()[flv]) / 1000 # convert to MeV
        ave_dict[time] = ave_

        # number of neutrinos emitted in this time bin, separated into energy bins
        egroup_ = [zero] # start with 0 neutrinos at 0 MeV bin
        for i in range_egroup:
            line = map(float, chunk[i].split())
            egroup_.append(line[flv-3])

            # Once, for the very first time bin, save the energy bins:
            if egroup_dict == {}:
                e_bins.append(line[1] / 1000) # energy of this bin (in MeV)

        egroup_dict[time] = egroup_

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
        ### Fill in spec_dict (energy spectrum) ###
        # Convert from (varying-size) energy bins to MeV^-1
        E_integ_ = 0
        spec_ = []
        egroup_ = egroup_dict[time] # list: number of neutrinos in different e_bins

        for (j, n) in enumerate(egroup_):
            if j == 0 or j == len(egroup_)-1:
                spec_.append(zero)
            else:
                spec_.append(n / (e_bins[j+1] - e_bins[j-1]))
                E_integ_ += (spec_[j-1] + spec_[j]) * (e_bins[j] - e_bins[j-1]) / 2

        # Normalise to 1
        spec_dict[time] = [x / E_integ_ for x in spec_]

        ### Fill in NL_dict (number luminosity) ###
        if i == 0:
            # TODO: how does this interact with data from wilson-nb.txt
            NL_dict[time] = zero
        else:
            prev_time = times[i-1]
            NL_dict[time] = (N_dict[time] - N_dict[prev_time]) / (time - prev_time)

        ### Fill in dNLde_dict (differential number luminosity) ###
        dNLde_ = [NL_dict[time] * spectrum for spectrum in spec_dict[time]]
        dNLde_dict[time] = dNLde_

        ### Get emission spectrum by log cubic spline interpolation ###
        log_group_e = [log10(e_bin) for e_bin in e_bins] # difference between base 10 and base e is just a constant factor of ln(10)
        log_dNLde_ = [log10(d) for d in dNLde_]
        interpolated_spec_dict[time] = interpolate.pchip(log_group_e, log_dNLde_)

    return None
