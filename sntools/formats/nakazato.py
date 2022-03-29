"""Parse Nakazato fluxes.

For simulations by Nakazato et al., Astrophys. J. Supp. 205 (2013) 2, arXiv:1210.6841
and Nakazato et al., Astrophys. J. 804 (2015) 75, arXiv:1503.01236.
Flux files are available at http://asphwww.ph.noda.tus.ac.jp/snn/index.html
"""

from scipy import interpolate
from sntools.formats import BaseFlux, get_endtime, get_raw_times, get_starttime


class Flux(BaseFlux):
    def parse_input(self, input, inflv, starttime, endtime):
        """Read simulations data from input file.

        Arguments:
        input -- prefix of file containing neutrino fluxes
        inflv -- neutrino flavor to consider
        starttime -- start time set by user via command line option (or None)
        endtime -- end time set by user via command line option (or None)
        """
        self.times = []
        self.dNLdE = {}

        with open(input) as infile:
            indata = [list(map(float, line.split())) for line in infile if not (line.startswith("#") or line.isspace())]

        # 21 lines per time bin
        chunks = [indata[21 * i: 21 * (i + 1)] for i in range(len(indata) // 21)]

        # input files contain information for e, eb & x in neighbouring columns,
        # so depending on the flavor, we might need an offset
        offset = {"e": 0, "eb": 1, "x": 2, "xb": 2}[inflv]

        # for each time bin, save data to dictionaries to look up later
        for chunk in chunks:
            # first line contains time
            time = chunk[0][0] * 1000  # convert to ms
            self.times.append(time)

            diff_number_flux, energy_mesh = [0], [0]  # flux = 0 at 0 MeV
            for bin_data in chunk[1:-1]:  # exclude first line (time) and last line (empty)
                number_flux = bin_data[2 + offset] / 1000.0  # convert 1/s to 1/ms
                luminosity = bin_data[5 + offset] * 624.151  # convert erg/s to MeV/ms
                diff_number_flux.append(number_flux)
                energy_mesh.append(luminosity / number_flux)

            self.dNLdE[time] = interpolate.pchip(energy_mesh, diff_number_flux)

        self.starttime = get_starttime(starttime, self.times[0])
        self.endtime = get_endtime(endtime, self.times[-1])
        self.raw_times = get_raw_times(self.times, self.starttime, self.endtime)

    def prepare_evt_gen(self, binned_t):
        """Pre-compute values necessary for event generation.

        Scipy/numpy are optimized for parallel operation on large arrays, making
        it orders of magnitude faster to pre-compute all values at one time
        instead of computing them lazily when needed.

        Argument:
        binned_t -- list of time bins for generating events
        """
        # unnecessary here; linear interpolation is fast enough to do it on demand
        return None

    def nu_emission(self, eNu, time):
        """Number of neutrinos emitted, as a function of energy.

        This is not yet the flux! The geometry factor 1/(4 pi r**2) is added later.
        Arguments:
        eNu -- neutrino energy
        time -- time ;)
        """
        # find previous/next time bin and perform linear interpolation
        for t_prev, t_next in zip(self.times[:-1], self.times[1:]):
            if time < t_next:
                break

        dNLdE_prev = self.dNLdE[t_prev](eNu)
        dNLdE_next = self.dNLdE[t_next](eNu)
        dNL = dNLdE_prev + (dNLdE_next - dNLdE_prev) * (time - t_prev) / (t_next - t_prev)

        return dNL
