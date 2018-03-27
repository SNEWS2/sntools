# sntools
Scripts for simulating a supernova neutrino burst in Hyper-Kamiokande.

### Required Input
Text file(s) containing information about neutrino fluxes. Multiple input formats are supported:

#### Garching format:
Three separate text files (one each for nu_e, anti-nu_e and nu_x - where nu_x stands for nu_mu or nu_tau or their respective antineutrinos), each containing (in this order) luminosity, mean energy and mean squared energy for different times - ideally in time steps of <1 ms – from a supernova simulation. See `sample-in.txt` for details.

#### Totani format:
Used by Totani et al. 1998, which is the baseline model in the Hyper-Kamiokande Design Report. This format is rather complex, see flux-totani.py for details.

#### Nakazato format:
Used by recent simulations of Nakazato et al., see http://asphwww.ph.noda.tus.ac.jp/snn/index.html
(Currently work in progress.)

### Interaction Channels
Currently the four main interaction channels in water Cherenkov detectors are supported:
inverse beta decay, elastic scattering on electrons and charged current interactions of nu_e and anti-nu_e on oxygen-16 nuclei.
For details, see the files in /interaction-channels/.

### Output:
A .kin file in the NUANCE format used by the /mygen/vecfile options in WCSim (see http://neutrino.phy.duke.edu/nuance-format/ for the format documentation).

### Typical Usage:
`python genevts.py --hierarchy=normal --channel=ibd -i infile -o outfile.kin`
This assumes the three input files (in Garching format) are named `infile_e.txt`, `infile_eb.txt` and `infile_x.txt`.

For a list of options, see `python genevts.py -h`.
