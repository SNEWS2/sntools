# sntools
Scripts for simulating a supernova neutrino burst in Hyper-Kamiokande.

Required input:
* information about neutrino fluxes (luminosity, mean energy and mean squared energy for different times - ideally in time steps of <1 ms) from a supernova simulation.
* comma-separated text file, one time step per line: `time (s), meanEnergy (MeV), meanSquaredEnergy (MeV^2), luminosity (erg/s)`

Output:
* .kin file in the NUANCE format used by the /mygen/vecfile options in WCSim
* for the format documentation, see http://neutrino.phy.duke.edu/nuance-format/

MakeKin_supernova.py was originally based on MakeKin.py, which is part of WCSim (https://github.com/WCSim/WCSim).
