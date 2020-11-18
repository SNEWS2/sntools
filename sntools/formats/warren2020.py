"""Parse Warren2020 fluxes.

Fluxes from https://zenodo.org/record/3952926 (DOI:10.5281/zenodo.3952926)
See https://arxiv.org/abs/1902.01340 and https://arxiv.org/abs/1912.03328
for description of the models.
"""

import h5py

from sntools.formats import gamma, get_starttime, get_endtime

flux = {}


def parse_input(input, inflv, starttime, endtime):
    """Read simulations data from input file.

    Arguments:
    input -- prefix of file containing neutrino fluxes
    inflv -- neutrino flavor to consider
    starttime -- start time set by user via command line option (or None)
    endtime -- end time set by user via command line option (or None)
    """

    f = h5py.File(input, 'r')
    for (t, r) in f['sim_data']['shock_radius']:
        if r > 1:
            tbounce = t * 1000  # convert to ms
            break

    starttime = get_starttime(starttime, 1000 * f['sim_data']['shock_radius'][0][0] - tbounce)
    endtime = get_endtime(endtime, 1000 * f['sim_data']['shock_radius'][-1][0] - tbounce)

    # Save flux data to dictionary to look up in nu_emission() below
    global flux
    flux = {}
    path = {'e': 'nue_data', 'eb': 'nuae_data', 'x': 'nux_data', 'xb': 'nux_data'}[inflv]
    for i, (t, lum) in enumerate(f[path]['lum']):
        t = 1000 * t - tbounce  # convert to time post-bounce in ms
        if (t < starttime - 30) or (t > endtime + 30):
            # Ignore data outside of the requested time span.
            continue

        lum *= 1e51 * 624.151  # convert from 10^51 erg/s to MeV/ms
        mean_e = f[path]['avg_energy'][i][1]
        mean_e_sq = f[path]['rms_energy'][i][1]**2

        flux[t] = (mean_e, mean_e_sq, lum)

    f.close()
    return (starttime, endtime, sorted(flux.keys()))


def prepare_evt_gen(binned_t):
    global flux
    gamma.flux = flux
    gamma.prepare_evt_gen(binned_t)
    flux = gamma.flux


def nu_emission(eNu, time):
    gamma.flux = flux
    return gamma.nu_emission(eNu, time)
