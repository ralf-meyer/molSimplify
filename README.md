![](./molSimplify/icons/logo.png)
[![Build Status](https://travis-ci.org/hjkgrp/molSimplify.svg?branch=master)](https://travis-ci.org/hjkgrp/molSimplify)

molSimplify is an open source toolkit for the automated, first-principles screening and discovery of new inorganic molecules and intermolecular complexes. molSimplify is developed by the [Kulik Group](http://hjkgrp.mit.edu) in the [Department of Chemical Engineering](http://web.mit.edu/cheme/) at [MIT](http://web.mit.edu). The software can generate a variety of coordination complexes of metals coordinated by ligands in a mono- or multi-dentate fashion. The code can build a coordination complex directly from a central atom or functionalize a more complex structure (e.g. a porphyrin or other metal-ligand complex) by including additional ligands or replacing existing ones. molSimplify also generates inter-molecular complexes for evaluating binding interactions and generating candidate reactants and intermediates for catalyst reaction mechanism screening.

## Installation

molSimplify is written in Python using the [Open Babel Toolkit](http://openbabel.org). The graphical user interface is written in PyQt5. We currently recommend installation via the [Conda](https://conda.io/docs/) package management system. Specific instructions are provided at the [Kulik group website](http://hjkgrp.mit.edu/content/new-installation-option-molsimplify). Additional installation options, including instructions for building the dependencies from source are provide in the [user guide](./molSimplify/Docs/Userguide/USER_GUIDE.pdf).

## Tutorials

A set of tutorials covering common use cases is available at the [Kulik group webpage](http://molsimplify.mit.edu).


## Citation [![DOI for Citing MDTraj](https://img.shields.io/badge/DOI-10.1002%2Fjcc.24437-blue.svg)](http://dx.doi.org/10.1002/jcc.24437)

molSimplify is research software. If you use it for work that results in a publication, please cite the following reference:

```
@article {molSimplify,
author = {Ioannidis, Efthymios I. and Gani, Terry Z. H. and Kulik, Heather J.},
title = {molSimplify: A toolkit for automating discovery in inorganic chemistry},
journal = {Journal of Computational Chemistry},
volume = {37},
number = {22},
issn = {1096-987X},
url = {http://dx.doi.org/10.1002/jcc.24437},
doi = {10.1002/jcc.24437},
pages = {2106--2117},
year = {2016},
}
```
