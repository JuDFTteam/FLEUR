Developing FLEUR
====================

The development effort for FLEUR is mainly hosted at [the Institute Quantum Theory of Materials @Forschungszentrum Juelich Germany](https://www.fz-juelich.de/pgi/pgi-1/EN).

GITLAB
------

The development process is performed using gitlab. You can access the [main gitlab page here](http://iffgit.fz-juelich.de/fleur/fleur).

If you checkout the code please be aware that there are several branches.

* The release branch contains the code of the last release published on the FLEUR webpage. You can not push to this branch directly. 
* You probably want to use the development branch to insert your changes. 
* If your changes are large, it might be a good idea to create your own branch first.

The changes you push to the gitlab will be tested by our CI directly: [![](https://iffgit.fz-juelich.de/fleur/fleur/badges/develop/pipeline.svg)](https://iffgit.fz-juelich.de/fleur/fleur/pipelines).

Doxygen
------
We use doxygen to create the documentation of the source. This can be [found here](https://fleur.iffgit.fz-juelich.de/fleur/doxygen).


Coverage
---------
The automatic tests of FLEUR cover only part of the source. [Here you find the analysis](https://fleur.iffgit.fz-juelich.de/fleur/coverage_html).


Further information
---------------
Some more information for developers are collected [here](developers.md).
