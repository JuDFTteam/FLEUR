#Differences between FLEUR v0.26 and v0.27

v0.27 is a major refactoring of FLEUR. Hence one can expect new errors, fixed bugs and in general slightly different results from this release. Here a list of factors that typically lead to changed results:

* When using inpgen new defaults have been set. This will usually result in more LOs and can change results.
* The behaviour of qpw_to_nmt is changed when using inp.xml. There is a new cutoff parameter for the l-expansion which is usually set to 0. This will change the core-tail correction and the generation of the starting charge density. Setting it to a large enough value (>=lmax) will switch back to the old behaviour.
* The previous precompiler option CPP_CORE is now standard. It can be switched in inp.xml.
* The precision of the eV-to-Hartree conversion factor increased. It was 27.2 in the old fleur version and now is 27.21138602. This factor is relevant for the generation of data files for DOS and bandstructure plots.
* The precision of the core solver has been increased. This will also have an influence on the starting density and the energy parameters.