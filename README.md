# sntools
Event generator for simulating supernova neutrino bursts in Hyper-Kamiokande (or other water Cherenkov detectors).


### Input
Text file(s) containing information about neutrino fluxes produced by the supernova.
sntools distinguishes between three flavours: nu_e, anti-nu_e and nu_x (where nu_x stands for nu_mu or nu_tau or their respective antineutrinos).
The following input formats are supported; see the source files in the `formats/` directory for details.

#### Garching format:
Text file containing time, mean energy, mean squared energy and luminosity. The flux is assumed to follow a Gamma distribution. See `sample-in.txt` for details.

#### Princeton format
Used in [recent simulations](https://arxiv.org/abs/1804.00689) by the Princeton group.

#### Nakazato format:
Used by recent simulations by Nakazato et al., [available for download here](http://asphwww.ph.noda.tus.ac.jp/snn/index.html).
#### Totani format
Used by Totani et al. 1998, which is the baseline model in the [Hyper-Kamiokande Design Report](https://arxiv.org/abs/1805.04163).

### Interaction Channels
Currently the four main interaction channels in water Cherenkov detectors are supported:
inverse beta decay, elastic scattering on electrons and charged current interactions of nu_e and anti-nu_e on oxygen-16 nuclei.
For details, see the files in `interaction_channels/`.


### Output
A .kin file in the NUANCE format used by the /mygen/vecfile options in WCSim. See [the format documentation](http://neutrino.phy.duke.edu/nuance-format/) for details.


### Typical Usage
```
python genevts.py sample-in.txt --format=garching -o outfile.kin --hierarchy=normal --channel=ibd
```

See
```
python genevts.py -h
```
for a full description of these and other options.

sntools currently uses Python 2.7, numpy 1.8 (or higher) and scipy 0.13 (or higher).
