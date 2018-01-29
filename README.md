# sntools
Scripts for simulating a supernova neutrino burst in Hyper-Kamiokande.

### Required input:
Three comma-separated text files (one each for nu_e, anti-nu_e and nu_x - where nu_x stands for nu_mu or nu_tau or their respective antineutrinos) containing information about neutrino fluxes (luminosity, mean energy and mean squared energy for different times - ideally in time steps of <1 ms) from a supernova simulation. See `sample-in.txt` for details.

### Output:
A .kin file in the NUANCE format used by the /mygen/vecfile options in WCSim (see http://neutrino.phy.duke.edu/nuance-format/ for the format documentation).

### Typical Usage:
`python genevts.py --hierarchy=normal --channel=ibd -i infile -o outfile.kin`
This assumes the three input files are named `infile_e.txt`, `infile_eb.txt` and `infile_x.txt`.

For a list of options, see `python genevts.py -h`.
