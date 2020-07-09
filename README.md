# sntools
sntools is a Monte Carlo event generator for supernova neutrino interactions.

Based on detailed time- and energy-dependent neutrino fluxes provided by various supernova models, it generates interactions within the detector volume and writes them to event files that can be used as an input for a full detector simulation.
sntools was originally developed for Hyper-Kamiokande and later extended to support different detectors and detector materials.

## Getting Started
First, make sure you have Python installed on your computer. (Either Python 2.7 or Python 3.x is fine.)

Then, in a terminal, run
```
pip install sntools
```
to install the latest version of sntools and all dependencies.
(sntools currently requires numpy 1.8 (or higher) and scipy 0.17 (or higher).)

You can then run
```
sntools -h
```
for a brief summary of all of sntools’ options.
A typical usage will look something like this:
```
sntools fluxes/intp2001.data --format nakazato --ordering normal --detector HyperK
```
This generates events for a supernova in Hyper-Kamiokande, assuming normal mass ordering and the neutrino flux given in the input file `fluxes/intp2001.data`, which is in the `nakazato` format.
(That sample file is included in the sntools repository. See [this web page](http://asphwww.ph.noda.tus.ac.jp/snn/index.html) for more information.)


## Input
Text file(s) containing information about neutrino fluxes produced by the supernova.
sntools distinguishes between three flavours: nu_e, anti-nu_e and nu_x (where nu_x stands for nu_mu or nu_tau or their respective antineutrinos).
The following input formats are supported; see the source files in the `formats/` directory for details.

#### Nakazato format
Used by recent simulations by the Nakazato group. Fluxes for 13 and 20 solar mass progenitors are included as `fluxes/intp1301.data` and `fluxes/intp2001.data`. [A description of the format and fluxes for more progenitors are available online.](http://asphwww.ph.noda.tus.ac.jp/snn/index.html)
If you use these included models in your work, please cite [Nakazato et al., ApJ Supp. 205 (2013) 2](https://arxiv.org/abs/1210.6841).

#### Gamma format
Text file containing time, mean energy, mean squared energy and luminosity. These parameters describe a Gamma distribution, [which is a good fit to the true spectrum](https://arxiv.org/abs/1211.3920). See `fluxes/sample-gamma.txt` for an unphysical sample file.

#### Princeton format
Used in [recent simulations](https://arxiv.org/abs/1804.00689) by the Princeton group.

#### Totani format
Used in [historical simulation by Totani et al.](https://arxiv.org/abs/astro-ph/9710203), which is the baseline model in the [Hyper-Kamiokande Design Report](https://arxiv.org/abs/1805.04163).


## Interaction Channels
sntools supports the main interaction channels in water and liquid scintillator.

For water Cherenkov detectors, like Hyper-Kamiokande, these are inverse beta decay, elastic scattering on electrons and charged-current interactions of nu_e and anti-nu_e on oxygen-16 nuclei.

For liquid scintillator detectors, these are inverse beta decay, elastic scattering on electrons, charged-current interactions of nu_e and anti-nu_e on carbon-12 nuclei and neutral-current interactions on carbon-12 nuclei.

Water-based liquid scintillator, a mixture of the two materials, is also supported.

For details, see the files in `interaction_channels/`.


## Output
A .kin file in the NUANCE format used by the /mygen/vecfile options in WCSim. See [the format documentation](http://neutrino.phy.duke.edu/nuance-format/) for details.
