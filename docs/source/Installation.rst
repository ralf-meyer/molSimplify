Installation
============

Regardless of install method, we recommend using anaconda to manage virtual environments. If you want a simpler install, you can download our pre-installed docker image, which has an out-of-the-box copy of molSimplify. We recommend installation via source for more advanced users. 


Installing from Source
----------------------

We recommend working in a virtual environment, even when installing from source. To get an up-to-date version of molSimplify, clone our repo:

::

    $ git clone https://github.com/hjkgrp/molSimplify.git

Go to the base folder for molSimplify, where you can create the conda environment from the yaml file (mols.yml). This step will help you get all the dependencies correct in a newly created conda environment named "mols_test". You can specify a different name for this environment at the first line of the yaml file (if you would like to change it from the default). Run the following commands:

::

    $ cd molSimplify/conda-envs
    $ conda env create -f mols.yml

This will take some time. After it is finished, activate your environment:

::

    $ conda activate mols_test

Now, go back to the root directory of molSimplify (where the setup.py file is located). Locally install molSimplify with pip. The period is necessary to specify the current directory.

::

    $ pip install -e .

To test your installation, you can run the command below at the root directory of molSimplify. You are good to go if all the tests are passed!

::

    $ python setup.py test


Using a Docker Image
--------------------

We also maintain an active docker image on [dockerhub](https://hub.docker.com/r/hjkgroup/molsimplify) for plug-and-play use.

For line by line instructions on an installation via docker, please visit [molSimplify installation webpage](http://hjkgrp.mit.edu/content/installing-molsimplify) of the Kulik group.