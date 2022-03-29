---
title: 'sntools: An event generator for supernova burst neutrinos'
tags:
  - Python
  - astrophysics
  - neutrino
  - supernova
authors:
  - name: Jost Migenda
    orcid: 0000-0002-5350-8049
    affiliation: "1,2" # (Multiple affiliations must be quoted)
  - name: Susan Cartwright # TODO: ORCID identifier
    affiliation: 2
  - name: Liz Kneale
    orcid: 0000-0002-4087-1244
    affiliation: 2
  - name: Matthew Malek # TODO: ORCID identifier
    affiliation: 2
  - name: Yan-Jie Schnellbach
    orcid: 0000-0002-6007-105X
    affiliation: 3
  - name: Owen Stone # TODO: ORCID identifier
    affiliation: 2
affiliations:
 - name: King’s College London
   index: 1
 - name: University of Sheffield
   index: 2
 - name: University of Liverpool
   index: 3
date: 14 April 2021
bibliography: documentation.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/1538-4357/abf7c4 # Update this with the DOI from AAS once you know it.
aas-journal: The Astrophysical Journal # The name of the AAS journal.

# This paper will be published in JOSS as DOI:10.21105/joss.2877
# See full review: https://github.com/openjournals/joss-reviews/issues/2877
---

# Summary
Supernovae are stellar explosions that produce many of the chemical elements
necessary for life to exist and form neutron stars and black holes.
However, despite millennia of observations and almost a century of in-depth
study, the explosion mechanism of supernovae is not yet well understood.

Observing the neutrino signal from the next galactic supernova can help solve
this puzzle.
To prepare for this once-in-a-lifetime event, it is essential to study how
various detectors respond to the neutrino fluxes from different supernovae.
sntools helps with this by providing a link between computer simulations of
supernovae and those of neutrino detectors:
From the neutrino fluxes predicted by a given supernova model, it generates
data sets of neutrino interactions, taking into account detailed,
energy-dependent cross-sections for all relevant interaction channels.
These data sets can then be used as input for a full detector simulation and
event reconstruction toolchain.

# Statement of Need
sntools is an event generator for supernova burst neutrinos which is written in
Python and makes extensive use of numpy [@Walt2011] and scipy [@Virtanen2020].
It currently supports multiple detector configurations using either water,
liquid scintillator or water-based liquid scintillator as a detection material.
It also supports several different input formats for neutrino fluxes from
various computer simulations.
New detector configurations, materials or input formats can be added easily.

sntools was initially developed to study supernova model discrimination
with Hyper-Kamiokande [@Migenda2019] and is also used to develop a supernova
DAQ system for the same experiment.
More recently, sntools was adapted by the WATCHMAN [@Askins2015] experiment
and for early studies for the THEIA detector concept [@Askins2020].

A few other software packages related to supernova neutrinos already exist.
SNOwGLoBES [@snowglobes] is widely used to compute event rates and energy
distributions for supernova burst neutrinos in various different detectors.
While it is an excellent tool for preliminary studies or quick comparisons of
different detector configurations, it uses simplified approximations for
detector effects like energy resolution or threshold.
It is an event rate calculator—not an event generator—and states in its own
documentation that it is “not intended to replace full detector simulations”.
In contrast, sntools is intended for use in conjunction with a full detector
simulation to perform more advanced, in-depth studies.

Software that performs a similar role to sntools was likely also created by
several different neutrino experiments.
However, it often is not publicly available or in use beyond the experiment it
was originally developed for.
One goal of sntools is to avoid such duplication of efforts in the future.
We welcome contributions, preferably by opening pull requests or issues on
the sntools GitHub repository.


# Acknowledgements
We want to thank members of the Hyper-Kamiokande software, astrophysics and
DAQ working groups for valuable discussions. We also want to thank Gabriel
Orebi Gann and Michael Wurm for advice on implementing the THEIA detector.


# References
