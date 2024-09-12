Installation instructions of FLEUR
=====================

To install FLEUR on your machine several options are possible:

1. Using the 'conda' package manager
2. Using the 'spack' package manager
3. Compiling from source

Using the conda package manager
-------------------------------
Using 'conda' can be an easy and painless way to get a working version of FLEUR on your machine. 
Please note:
* The FLEUR version you obtain is not necessary the newest and
* ***most importantly*** the executable will not be optimized for your system and will not be HPC ready.
* we only recommend this path if you want to try out FLEUR without much hassle but strongly recommend a proper install as described below for production runs.

If you have a 'conda' on your system installing FLEUR should work with the command:

`conda install conda-forge::fleur`

Using the spack package manager
-------------------------------
The 'spack' package manager can offer a very efficient path to install a version of FLEUR optimized to your system including many dependencies you might want to use.

In order to install FLEUR using 'spack' you should first:
* install 'spack' as described here https://spack.readthedocs.io/en/latest/getting_started.html .
* clone the FLEUR-sources from the git repo: `git clone https://iffgit.fz-juelich.de/fleur/fleur`

Now you might create a spack environment to install FLEUR. There are some examples in the `packaging/spack/environments` directories which you might want to use or use as starting points.
E.g. you might want to start by `spack env create FLEUR-INSTALL packaging/spack/environments/oneapi-hdf-elsi.yaml`. 
*** Please note that we have some modified spack recipies in `packaging/spack/repo` so you should use this repo in you environment. For the examples provided this is done if you set the environment variable FLEUR_SOURCE to the directory in which you cloned FLEUR.***

Then activate and adjust your environment and install FLEUR as described in the spack users guide: https://spack.readthedocs.io/en/latest/environments.html

Compiling from source
--------------

Currently the standard path to install FLEUR is by installing it from the source. Please visit https://www.flapw.de/MaX-7.0/documentation/installation/ for information.
