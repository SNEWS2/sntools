# sntools
sntools is a Monte Carlo event generator for supernova neutrino interactions.

It is a command line application that uses detailed time- and energy-dependent neutrino fluxes (provided by various supernova models) to generate interactions within the detector volume and write them to event files that can be used as an input for a full detector simulation.
sntools was originally developed for Hyper-Kamiokande and later extended to support different detectors and detector materials.

This README file should give a brief overview over sntools and help you get started. For more information, see the [full documentation for each release on GitHub](https://github.com/JostMigenda/sntools/releases).

## Getting Started
### Installation Instructions
First, make sure you have Python installed on your computer. (Either Python 2.7 or Python 3.x is fine.)

Then, in a terminal, run
```
pip install sntools
```
to install the latest version of sntools and all dependencies.
(sntools currently requires at least numpy 1.8, scipy 0.17 and h5py 2.10, but newer versions are recommended.)

Finally, run
```
python -c 'import sntools; sntools.setup()'
```
to check whether sntools is working and to download the sample input file used below.

### Example Usage

To generate your first supernova neutrino events with sntools, open a terminal window and type in the following command:
```
sntools fluxes/intp2001.data --format nakazato --endtime 50
```
This generates events for the first 50 ms of a supernova based on the neutrino flux given in the input file `fluxes/intp2001.data` (which is in the `nakazato` format) and writes them to a file named `outfile.kin` in the current directory.
(This sample input file is included in the sntools repository. See [this web page](http://asphwww.ph.noda.tus.ac.jp/snn/index.html) for more information.)

A more realistic usage that demonstrates more of sntools’ capabilities looks like this:
```
sntools fluxes/intp2001.data --format nakazato --channel es --detector SuperK --distance 20 --verbose --output intp2001es.kin
```
This uses the neutrino flux given in the same input file to generate neutrino-electron elastic scattering events in Super-Kamiokande for a supernova at 20 kpc distance. It also produces more verbose output, which lets you see that sntools generates events separately for different neutrino flavours (which have different fluxes and cross sections), before merging them into an output file named `intp2001es.kin`.

You can also run
```
sntools -h
```
to get an overview over all of sntools’ options.


## Input
Text file(s) containing information about neutrino fluxes produced by the supernova.
sntools distinguishes between three flavours: nu_e, anti-nu_e and nu_x (where nu_x stands for nu_mu or nu_tau or their respective antineutrinos).
The following input formats are supported.

#### Nakazato format
Used by recent simulations by the Nakazato group. Fluxes for 13 and 20 solar mass progenitors are included as `fluxes/intp1301.data` and `fluxes/intp2001.data`. [A description of the format and fluxes for more progenitors are available online.](http://asphwww.ph.noda.tus.ac.jp/snn/index.html)
If you use these included models in your work, please cite [Nakazato et al., ApJ Supp. 205 (2013) 2](https://arxiv.org/abs/1210.6841).

#### Gamma format
Text file containing time, mean energy, mean squared energy and luminosity. These parameters describe a Gamma distribution, [which is a good fit to the true spectrum](https://arxiv.org/abs/1211.3920). See `fluxes/sample-gamma.txt` for an unphysical sample file.

#### Warren2020 format
Similar to the Gamma format, but in an HDF5 file instead of plain text. Data [available online](https://zenodo.org/record/3952926).

#### Princeton format
Used in [recent simulations](https://arxiv.org/abs/1804.00689) by the Princeton group.

#### Totani format
Used in [historical simulation by Totani et al.](https://arxiv.org/abs/astro-ph/9710203), which is the baseline model in the [Hyper-Kamiokande Design Report](https://arxiv.org/abs/1805.04163).


## Interaction Channels
sntools supports the main interaction channels in water and liquid scintillator.

For water Cherenkov detectors, like Hyper-Kamiokande, these are inverse beta decay, elastic scattering on electrons and charged-current interactions of nu_e and anti-nu_e on oxygen-16 nuclei.

For liquid scintillator detectors, these are inverse beta decay, elastic scattering on electrons, charged-current interactions of nu_e and anti-nu_e on carbon-12 nuclei and neutral-current interactions on carbon-12 nuclei.

Water-based liquid scintillator, a mixture of the two materials, is also supported.


## Output
A text file in the NUANCE format (used by the `/mygen/vecfile` options in WCSim) or the RATPAC format.

## Support and Contributing
To report problems or ask for help, [open an issue on GitHub](https://github.com/JostMigenda/sntools/issues) or email the [lead developer](https://github.com/JostMigenda/).

Contributions to sntools are welcome! See instructions in the [full documentation](https://github.com/JostMigenda/sntools/releases) for help on common ways to extend sntools (e. g. by adding a new detector, input format or interaction channel) and please submit a pull request with your contributions!

Please note that this project is released with a Contributor Code of Conduct (see `CODE_OF_CONDUCT.md`). By participating in this project you agree to abide by its terms.
