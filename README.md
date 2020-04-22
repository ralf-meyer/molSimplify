![](./molSimplify/icons/logo.png)
[![Build Status](https://travis-ci.org/hjkgrp/molSimplify.svg?branch=master)](https://travis-ci.org/hjkgrp/molSimplify)

molSimplify is an open source toolkit for the automated, first-principles screening and discovery of new inorganic molecules and intermolecular complexes. molSimplify is developed by the [Kulik Group](http://hjkgrp.mit.edu) in the [Department of Chemical Engineering](http://web.mit.edu/cheme/) at [MIT](http://web.mit.edu). The software can generate a variety of coordination complexes of metals coordinated by ligands in a mono- or multi-dentate fashion. The code can build a coordination complex directly from a central atom or functionalize a more complex structure (e.g. a porphyrin or other metal-ligand complex) by including additional ligands or replacing existing ones. molSimplify also generates inter-molecular complexes for evaluating binding interactions and generating candidate reactants and intermediates for catalyst reaction mechanism screening. molSimplify also ships neural network models that can predict the [metal-ligand bond lengths](https://pubs.rsc.org/en/content/articlehtml/2017/sc/c7sc01247k), [spin-splitting energy](https://pubs.acs.org/doi/abs/10.1021/acs.jpca.7b08750), [frontier orbital energies](https://pubs.acs.org/doi/abs/10.1021/acs.iecr.8b04015), and [simulation outcomes](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.9b00057) for octahedral transition metal complexes. See the Tutorials at the [Kulik group webpage](http://hjkgrp.mit.edu/molSimplify-tutorials) for a more complete list of jobs molSimplify can do.

## Installation

### via conda
We currently recommend installation via the [Conda](https://conda.io/docs/) package management system.
1. Prerequisite: have [Anaconda or miniconda](https://www.anaconda.com/distribution/) installed on your system.

2. Clone molSimplify source from github.

   ```bash
   git clone https://github.com/hjkgrp/molSimplify.git
   ```
   
3. Go to the folder root folder for molSimplify, create the conda environment from the yaml file (mols.yml). This step will help you get all the dependecies correct in a newly created conda environment named "mols_test". You can specify a different name for this environment at the first line of the yaml file.

   ```bash
   cd <root directory to molSimplify>/conda-envs
   conda env create -f mols.yml
   ```
4. Activate the conda environment you just created. Go back to the root directory of molSimplify (where the setup.py file locates). Local install with pip.
   ```bash
   cd <root directory to molSimplify>
   pip install -e .
   ```
5. To test your installation, you can run the command below at the root directory of molSimplify. You are good to go if all the tests are passed!
   ```bash
   python setup.py test
   ```

### via docker
We also maintain an active [docker image on dockerhub](https://hub.docker.com/repository/docker/hjkgroup/molsimplify) for plug-and-play use. 

For line by line instructions on an installation via docker, please visit [molSimplify installation webpage of Kulik group](http://hjkgrp.mit.edu/content/installing-molsimplify).

## Tutorials

A set of tutorials covering common use cases is available at the [Kulik group webpage](http://hjkgrp.mit.edu/molSimplify-tutorials).


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

