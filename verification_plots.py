#!usr/bin/python
"""
Plot energy and event rate of generated events.

This is useful to check results of event generation phase. See
`python verification_plots.py -h` for usage information.
"""

import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser

parser = OptionParser()

optdefault = "tmp_ibd_eb.txt"
parser.add_option("-i", "--input", dest="input",
                  help="File containing events to plot. Default: '%s'." \
                      % optdefault,
                  metavar="FILENAME",
                  default=optdefault)

parser.add_option("-a", "--all", dest="show_all",
                  help="Show all plots. Off by default.",
                  default=False, action="store_true")

parser.add_option("--spectrum", dest="show_spectrum",
                  help="Show time-integrated energy spectrum. Off by default.",
                  default=False, action="store_true")

parser.add_option("--meanene", dest="show_meanene",
                  help="Show mean energy in each time bin. Off by default.",
                  default=False, action="store_true")

parser.add_option("--evtcount", dest="show_evtcount",
                  help="Show event count in each time bin. Off by default.",
                  default=False, action="store_true")

(options, args) = parser.parse_args()

if options.show_all:
    options.show_spectrum = True
    options.show_meanene  = True
    options.show_evtcount = True

# Get energy and time of generated events.
energy = []
time = []
with open(options.input) as data: # for elastic scattering: "tmp_es_ee.txt"
    for line in data:
        t, _, e, _, _, _ = line.split(",")
        energy.append(float(e))
        time.append(float(t))

# Get event counts and mean energies for each time bin.
bin_width = 10 # time interval in ms
t_max = 320
bin_nr = t_max / bin_width
total_energy = [0] * bin_nr
n_evt = [0] * bin_nr

for (t, energy) in zip(time, energy):
    if t >= t_max or t < 0:
        continue

    j = int(t / bin_width)
    total_energy[j] += + energy
    n_evt[j] += 1

mean_energy = []
for j in range(bin_nr):
    if n_evt[j] == 0:
        mean_energy.append(0)
    else:
        mean_energy.append(total_energy[j] / n_evt[j])


# import data from pre-processed supernova simulation data (in Garching format)
# sim_mean_energies = []
# sim_times = []
# with open("infile_eb.txt") as data2:
#     for line in data2:
#         t, mean_energy, _, _ = line.split(",")
#         sim_times.append(float(t) * 1000) # convert to ms
#         sim_mean_energies.append(float(mean_energy))

if options.show_spectrum:
    # plot energy spectrum
    plt.title("Energy spectrum")
    plt.ylabel("Events")
    plt.xlabel("Energy (MeV)")
    plt.hist(energy, bins=50, histtype="step")
    plt.show()

    # as above, but with logarithmic y-axis
    plt.title("Energy spectrum (log)")
    plt.ylabel("Events")
    plt.xlabel("Energy (MeV)")
    plt.hist(energy, bins=50, histtype="step", log=(True))
    plt.show()

if options.show_meanene:
    # plot energy against time for charged particles (and neutrinos from supernova simulation)
    plt.title("Mean energy of charged particle")
    plt.ylabel("Mean energy (MeV)")
    plt.xlabel("Time (ms)")
    plt.plot(range(0, t_max, bin_width), mean_energy)
    # plt.plot(sim_times, sim_mean_energies, 'r') # adds a plot of the neutrino energies direct from the supernova simulation data
    plt.show()

if options.show_evtcount:
    plt.title("Event count per bin")
    plt.ylabel("Events")
    plt.xlabel("Time (ms)")
    plt.plot(range(0, t_max, bin_width), n_evt)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
