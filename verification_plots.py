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
energies = []
times = []
with open(options.input) as data: # for elastic scattering: "tmp_es_ee.txt"
    for line in data:
        t, _, energy, _, _, _ = line.split(",")
        energies.append(float(energy))
        times.append(float(t))

# Get event counts and mean energies for each time bin.
bin_width = 2 # time interval in ms
bin_nr = np.arange(1, 500/bin_width, 1)
mean_energies = []
n_evts = []

for i in bin_nr:
    total_energy = 0
    n_evt = 0
    time = 15 + (i*bin_width)
    boundsMin = time - bin_width
    boundsMax = time

    for (i, time) in enumerate(times):
        if boundsMin <= time < boundsMax:
            total_energy = total_energy + energies[i]
            n_evt = n_evt + 1

    if n_evt == 0:
        mean_energies.append(total_energy)
    else:
        mean_energies.append(total_energy/n_evt)
    n_evts.append(n_evt) # create a list of number of events in each time bin of specified width

# import data from pre-processed supernova simulation data
sim_mean_energies = []
sim_times = []
with open("infile_eb.txt") as data2: # for elastic scattering: "infile_e.txt"
    for line in data2:
        t, mean_energy, _, _ = line.split(",")
        sim_times.append(float(t) * 1000) # convert to ms
        sim_mean_energies.append(float(mean_energy))

if options.show_spectrum:
    # plot energy spectrum
    plt.title("Energy spectrum")
    plt.ylabel("Events")
    plt.xlabel("Energy (MeV)")
    plt.hist(energies, bins=50, histtype="step")
    plt.show()

    # as above, but with logarithmic y-axis
    plt.title("Energy spectrum (log)")
    plt.ylabel("Events")
    plt.xlabel("Energy (MeV)")
    plt.hist(energies, bins=50, histtype="step", log=(True))
    plt.show()

if options.show_meanene:
    # plot energy against time for charged particles (and neutrinos from supernova simulation)
    plt.title("Mean energy of charged particle")
    plt.ylabel("Mean energy (MeV)")
    plt.xlabel("Time (ms)")
    plt.plot(bin_width*bin_nr+15, mean_energies)
    # plt.plot(sim_times, sim_mean_energies, 'r') # adds a plot of the neutrino energies direct from the supernova simulation data
    plt.show()

if options.show_evtcount:
    plt.title("Event count per bin")
    plt.ylabel("Events")
    plt.xlabel("Time (ms)")
    plt.plot(bin_width*bin_nr+15, n_evts)
    plt.show()
