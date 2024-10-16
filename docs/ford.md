exclude_dir: ../docs
exclude_dir: ../cmake
exclude_dir: ../tests
exclude_dir: ../build
exclude_dir: ../build.debug
exclude_dir: ../external
exclude: jpSetupDynMatDeprecated_mod.F90
exclude: soc_or_ssdw.f90
exclude: jpSternhHF_mod.F90
include: ../include/
include: ../init/compileinfo.h
src_dir: ../
output_dir: ./ford-doc
project: FLEUR
project_website: https://www.flapw.de
project_gitlab: https://iffgit.fz-juelich.de/fleur/fleur
project_url: https://fleur.iffgit.fz-juelich.de/fleur/ford
summary: The FLEUR code: All-electron full-potential augmented plane-wave method for DFT.
author: FLEUR development team
author_description: FLEUR is developed by an active community mostly based at the [Forschungszentrum Jülich](http://www.fz-juelich.de/pgi/pgi-1/EN/Home/home_node.html).
email: d.wortmann@fz-juelich.de
predocmark: >>
docmark_alt: !#!
media_dir: ./media
pages_dir: ./ford-pages
display: public
         protected
         private
source: true
graph: true
graph_maxnodes: 100
graph_maxdepth: 10
search: true
preprocess: false
macro: TEST
       LOGIC=.true.
license: mit

This is the automatic source code generation using [FORD](https://github.com/Fortran-FOSS-Programmers/ford). Please use it to explore the code. A good start point might be the "main fleur routine" in the module [[m_fleur]].

Further useful links:

* The [Homepage](https://www.flapw.de)
* The [Git-Repo](https://iffgit.fz-juelch.de/fleur/fleur)
