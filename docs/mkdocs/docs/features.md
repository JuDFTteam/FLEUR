Description of the FLEUR codes
================
The FLEUR code family is a program package for calculating ground-state as well as excited-state properties of solids. It is based on the full-potential linearized augmented-plane-wave (FLAPW) method [1-4]. The strength of the FLEUR code [5,6] lies in applications to bulk, semi-infinite, two- and one-dimensional solids [7], solids of all chemical elements of the periodic table, solids with complex open structures, low symmetry, with complex non-collinear magnetism [8] in combination with spin-orbit interaction [9,10], external electric fields, and the treatment of spin-dependent transport properties [11,12]. It is an all-electron method and thus treats core and valence electrons and can deal with hyperfine properties. The inclusion of local orbitals allows a systematic extension of the LAPW basis that enables a precise treatment of semicore states [13], unoccupied states [14,15], and an elimination of the linearization error in general [16]. A large variety of local and semi-local (GGA) exchange and correlation functionals are implemented, including the LDA+U approach. In recent years the code has been developed further to make contact to electronically complex materials. Hybrid functionals [17,18] and the optimized-effective-potential (OEP) method [15,19] have been implemented. Wannier functions [20] can be constructed to make contact to realistic model Hamiltonians. Excitations can be treated on the basis of the GW approximation [21,22] and ladder diagrams are included to compute spin-wave excitations [23]. The Hubbard U can be calculated in the constrained random phase approximation (cRPA) [24].

Literature:

1. O.K. Andersen, "Linear methods in band theory", [ Phys. Rev. B 12, 3060 (1975)](http://dx.doi.org/10.1103/PhysRevB.12.3060 )
2. D. D. Koelling and G. O. Arbman, Use of energy derivative of the radial solution in an augmented plane wave method: application to copper, [J. Phys. F: Metal Phys. 5, 2041 (1975)](http://dx.doi.org/10.1088/0305-4608/5/11/016 )
3. E. Wimmer,A.J. Freeman, H. Krakauer, and M. Weinert, "Full-potential self-consistent linearized-augmented-plane-wave method for calculating the electronic structure of molecules and surfaces: O2 molecule", [ Phys. Rev. B 24, 864 (1981)](http://dx.doi.org/10.1103/PhysRevB.24.864 )
4. M. Weinert, E. Wimmer, and A.J. Freeman, Total-energy all-electron density functional method for bulk solids and surfaces, [ Phys. Rev. B 26, 4571 (1982)](http://dx.doi.org/10.1103/PhysRevB.26.4571 )
5. S. Blügel and G. Bihlmayer, "[ Full-Potential Linearized Augmented Planewave Method](http://www2.fz-juelich.de/nic-series/volume31/bluegel.pdf )", in Computational Nanoscience: Do It Yourself! edited by J. Grotendorst, S. Blügel, and D. Marx, NIC Series Vol. 31, p. 85 (John von Neumann Institute for Computing, Jülich, 2006)
6. http://www.flapw.de
7. Y. Mokrousov, G. Bihlmayer, and S. Blügel, "A full-potential linearized augmented planewave method for one-dimensional systems: gold nanowire and iron monowires in a gold tube", [ Phys. Rev. B. 72, 045402 (2005)](http://dx.doi.org/10.1103/PhysRevB.72.045402 )
8. Ph. Kurz, F. Foerster, L.Nordström, G. Bihlmayer, and S. Blügel, "Ab initio treatment of non-collinear magnets with the full-potential linearized augmente planewave method", [ Phys. Rev. B 69, 024415 (2004)](http://dx.doi.org/10.1103/PhysRevB.69.024415 )
9. M. Heide, G. Bihlmayer, and S. Blügel, "Describing Dzyaloshinskii-Moriya spirals from first principles", [ Physica B 404, 2678 (2009)](http://dx.doi.org/10.1016/j.physb.2009.06.070 )
10. B. Zimmermann, M. Heide, G. Bihlmayer, and S. Blügel, "First-principles analysis of a homochiral cycloidal magnetic structure in a monolayer Cr on W(110)", [ Phys. Rev. B 90, 115427 (2014)](http://dx.doi.org/10.1103/PhysRevB.90.115427 )
11. D. Wortmann, H. Ishida, and S. Blügel, "Ab initio Green-function formulation of the transfer matrix: Application to complex bandstructures", [ Phys. Rev. B 65, 165103 (2002)](http://dx.doi.org/10.1103/PhysRevB.65.165103 )
12. D. Wortmann, H. Ishida, and S. Blügel, "Embedded Green-function approach to the ballistic electron transport through an interface", [ Phys. Rev. B 66, 075113 (2002)](http://dx.doi.org/10.1103/PhysRevB.66.075113 )
13. D. Singh, "Ground-state properties of lanthanum: Treatment of extended-core states", [ Phys. Rev. B 43, 6388 (1991)](http://dx.doi.org/10.1103/PhysRevB.43.6388 )
14. C. Friedrich, A. Schindlmayr, S, Blügel, and T. Kotani, "Elimination of the linearization error in GW calculations based on the linearized augmented-plane-wave method", [ Phys. Rev. B 74, 045104 (2006)](http://dx.doi.org/10.1103/PhysRevB.74.045104 )
15. M. Betzinger, C. Friedrich, S. Blügel, and A. Görling, "Local exact exchange potentials within the all-electron FLAPW method and a comparison with pseudopotential results", [ Phys. Rev. B 83, 045105 (2011)](http://dx.doi.org/10.1103/PhysRevB.83.045105 )
16. G. Michalicek, M. Betzinger, C. Friedrich, and S. Blügel, "Elimination of the linearization error and improved basis-set convergence within the FLAPW method", [ Comp. Phys. Commun. 184, 2670 (2013)](http://dx.doi.org/10.1016/j.cpc.2013.07.002 )
17. M. Betzinger, C. Friedrich, and S. Blügel, "Hybrid functionals within the all-electron FLAPW method: implementation and applications of PBE0", [ Phys. Rev. B 81, 195117 (2010)](http://dx.doi.org/10.1103/PhysRevB.81.195117 )
18. M. Schlipf, M. Betzinger, C. Friedrich, M. Ležai&#263;, and S. Blügel, "HSE hybrid functional within the FLAPW method and its application to GdN", [ Phys. Rev. B 84, 125142 (2011)](http://dx.doi.org/10.1103/PhysRevB.84.125142 )
19. M. Betzinger, C. Friedrich, A. Görling, and S. Blügel, "Precise response functions in all-electron methods: Application to the optimized-effective-potential approach", [ Phys. Rev. B 85, 245124 (2012)](http://dx.doi.org/10.1103/PhysRevB.85.245124 )
20. F. Freimuth, Y. Mokrousov, D. Wortmann, S. Heinze, and S. Blügel, "Maximally Localized Wannier Functions within the FLAPW formalism", [ Phys. Rev. B. 78, 035120 (2008)](http://dx.doi.org/10.1103/PhysRevB.78.035120 )
21. C. Friedrich, S. Blügel, and A. Schindlmayr, "Efficient implementation of the GW approximation within the all-electron FLAPW method", [ Phys. Rev. B 81, 125102 (2010)](http://dx.doi.org/10.1103/PhysRevB.81.125102 )
22. C. Friedrich, S. Blügel, and A. Schindlmayr, "Efficient calculation of the Coulomb matrix and its expansion around k=0 within the FLAPW method", [ Comp. Phys. Comm. 180, 347 (2009)](http://dx.doi.org/10.1016/j.cpc.2008.10.009 )
23. E. &#350;a&#351;&#305;o&#287;lu, A. Schindlmayr, Ch. Friedrich, F. Freimuth, and S. Blügel, "Wannier-function approach to spin excitations in solids", [ Phys. Rev. B 81, 054434 (2010)](http://dx.doi.org/10.1103/PhysRevB.81.054434 )
24. E. &#350;a&#351;&#305;o&#287;lu, C. Friedrich, and S. Blügel, "Effective Coulomb interaction in transition metals from constrained random-phase approximation", [ Phys. Rev. B 83, 121101(R) (2011)](http://dx.doi.org/10.1103/PhysRevB.83.121101 )

Features of FLEUR
==============

The FLEUR code allows you to investigate structural, electronic and magnetic 
properties of periodic systems, in bulk (3D), film (2D) and wire (1D) 
geometry. Furthermore, it provides the necessary input for the calculation 
of non-periodic systems (semi-infinite crystals or transport geometries) within 
the G-Fleur code, or for the calculation of excited state properties.

FLEUR is based on density functional theory ([DFT](http://www2.fz-juelich.de/nic-series/volume31/jones.pdf)) and is an implementation of the 
full-potential linearized augmented planewave ([FLAPW](http://www2.fz-juelich.de/nic-series/volume31/bluegel.pdf)) method, which 

* is a highly-precise all electron method
* has a basis set equally suited for open and close-packed systems
* is suitable for elements from the whole periodic table

Among other things, FLEUR allows to calculate

* structural and magnetic ground state properties
* electronic properties like bandstructures, densities of states etc.
* charge densities, field gradients, or hyperfine fields

Although FLEUR calculations can be performed for all kinds of materials, it 
is especially suited for:

* magnetic systems (collinear or non-collinear)
* open systems (surfaces, wires, nanostructures)
* transition metals, lanthanides, actinides
