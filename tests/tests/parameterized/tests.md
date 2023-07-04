Regression tests for FLEUR
================
This files contains the regression tests. A test definition consists of:

* A short description
* The directory of the input files. (These have to be sorted into the subdirectories named according to the testset.)
* The pytest markers for the tests

Only tests with a "+" in the first column are actually executed.

To create a test please append the corrsponding table.

Testset: Basic
------

||Description|directory name|marks|Remarks
|-|-----------|--------------|------|----
|+|Bulk Cu one iteration,bandstructure|basic/CuBand|
|+|Bulk Cu one iteration|basic/CuBulk|fast,bulk
|+|Bulk Cu one iteration,DOS|basic/CuDOS|fast,bulk,dos
||Bulk Cu one iteration,DOS,Orbital decomp.|basic/CuOrb|fast,bulk,dos,orbcomp
||Bulk Co, DOS,MCD|basic/CoMCD|bulk,dos,mcd
|+|Bulk Co, bandstructure, unfolding|basic/CoUnfold|band,bulk
|+|Bulk Fe, Kerker preconditioner|basic/Fe_Kerker|bulk
||Bulk Fe fcc with relativistic core solver|basic/Fe_fcc_kcrel|
|+|Si with LOs|basic/SiLO|bulk
|+|Bulk PTO|basic/PTO|bulk
|+|Bulk PTO, SOC|basic/PTO-SOC|bulk,soc
|+|Bulk Fe, Tetrahedon method|basic/Fe_Tetra_noSYM|bulk
||LDA+U with Around mean field double counting and magnetism|basic/NiOldaUAMF|bulk,ldau
|+|LDA+U with fully localized double counting and magnetism|basic/NiOldaUFLL|bulk,ldau
|+|Crystal field output|basic/CrystalFieldOutput|bulk

Testset: Films
----------------
||Description|directory name|marks|Remarks
|-|-----------|--------------|------|----
|+|Fe Monolayer SOC|film/Fe_1l_SOC|soc
|+|Fe Monolayer|film/Fe_1l|
|+|Fe Monolayer Triangular method|film/Fe_1l_Tria|
|+|Si Film, plotting|film/SiFilmPlot|
|+|Si Film, plotting,slicing|film/SiFilmSlicePlot|



Testset: Forces
------------
||Description|directory name|marks|Remarks
|-|-----------|--------------|------|---
|+|Bulk GaAs, Relaxation, LDA+U|forces/GaAsMultiUForce|bulk,ldau
|+|Bulk VO2, Relaxation         |forces/VO2_forces|bulk
|+|Bulk VO2, Relaxation, different force levels         |forces/VO2_force_levels|bulk
|+|Bulk H2O, Relaxtion using BFGS|forces/H2ORelaxBFGS|bulk

Testset: Noco
----------
||Description|directory name|marks|Remarks
|-|-----------|--------------|------|-----
|+|Fe bct,noco, SOCL|noco/Fe_bct_SOC|bulk,soc
|+|Fe bct,noco, LOs|noco/Fe_bct_LO|bulk
|+|Fe bct,noco|noco/Fe_bct|bulk
|+|Fe bct,noco,libxc|noco/Fe_bct_LibXC|bulk,libxc
||Noco, one atom, mag. in x direction|noco/1atx|bulk
||Noco, one atom, mag. in y direction|noco/1aty|bulk
||Noco, one atom, mag. in z direction|noco/1atz|bulk
||Noco, one atom, mag. in non-sym direction|noco/1at|bulk
||Noco,SOC, one atom, mag. in x direction|noco/1atSOCx|bulk,soc
||Noco,SOC, one atom, mag. in y direction|noco/1atSOCy|bulk,soc
||Noco,SOC, one atom, mag. in z direction|noco/1atSOCz|bulk,soc
||Noco,SOC, one atom, mag. in non-sym direction|noco/1atSOC|bulk,soc
|+|Fe fcc spin-spiral|noco/Fe_fcc|bulk
|+|Iron LO's and SOC test in FFN|noco/FeFFNLOsSOC|bulk,soc,hdf
|+|Fe monolayer fcc (110): spin spiral energy|noco/FePt_film_SSFT|film,spinspiral,forcetheorem
|+|Fe monolayer fcc (110): spin spiral energy, with LO|noco/FePt_film_SSFT_LO|film,spinspiral,forcetheorem
||Fe bcc, Flipcdn and noco in MT,x-dir|noco/Fe_bcc_FlipcdnXLDA|bulk|produces warnings
||Fe bcc, Flipcdn and noco in MT,y-dir|noco/Fe_bcc_FlipcdnYLDA|bulk|produces warnings
||relaxation feature of FFN in the MT|noco/RelaxMT|bulk,hdf

Testset: Experimental
----------
||Description|directory name|marks|Remarks
|-|-----------|--------------|------|-----
|+|Test of the orbital polarization correction|extra/Fe_bcc_OPC|bulk,soc
|+|Sourcefree magnetism and magnetization scaling|extra/Fe_bcc_SF_LDA|bulk,hdf
|+|Bulk Cu one iteration|extra/CuBulkLibXC|libxc,bulk
||Bulk Al one iteration, LibXC|extra/Al_libxc_PBE|bulk,libxc
||Test of GW interface 1 |extra/gw1Interface|bulk|inp.xml files too old
||Test of GW interface 2 |extra/gw2Interface|bulk|inp.xml files too old
|+|Sm jDOS decomposition|extra/SmAtomjDOS|bulk,dos
||C: simple test for the Wannier code|extra/Cwann|bulk,wannier
||TiO2 EELS spectrum|extra/TiO2eels|bulk,eels|inp.xml too old
|+|Hubbard1 using SOC|extra/Gd_Hubbard1|bulk,edsolver
|+|Hubbard1 without sym|extra/Gd_Hubbard1_noSYM|bulk,edsolver
||diamond for one k-point with scan|extra/Diamond_SCAN|bulk,libxc
|+|3D vector plots of the magnetization|extra/PlotOnlyMT|bulk,plot,hdf
|+|density and potential plots, vector plots|extra/PlotDenandPot|bulk,plot,hdf


Testset: Hybrid-Functionals
----------
||Description|directory name|marks|Remarks
|-|-----------|--------------|------|---
||GaAs PBE0|hybrid/GaAsHybridPBE0|bulk|runs too long
||Fe PBE0|hybrid/FeHybridPBE0|bulk|runs too long
||KCl PBE0|hybrid/KClHybridPBE0|bulk|runs too long
||Mn Noinversion, PBE0|hybrid/MnHybridNoinv|bulk|runs too long



Testset: Greenfunctions
----------
||Description|directory name|marks|Remarks
|-|-----------|--------------|------|---
|+|Fe bcc Green's function|greens/Fe_bcc_GreensFunction|bulk
|+|Fe Monolayer Green's function|greens/Fe_1l_GreensFunction|film
|+|Greens Function MultiContour|greens/GreensFunction_MultiContour|bulk
|+|Fe bcc Green's function Radial|greens/GreensFunctionRadial|bulk
||Fe bcc Green's function Radial with local orbitals|greens/GreensFunctionRadial_LO|bulk
|+|Ho atom Green's function|greens/GreensFunction_HoAtom_SQA_theta|bulk
|+|Ho atom Green's function|greens/GreensFunction_HoAtom_SQA_phi|bulk
|+|Ho atom Green's function|greens/GreensFunction_rotated_SQA_noco|bulk
|+|Fe bcc Green's function Radial Noco spin offdiagonal|greens/GreensFunction_mperp_xdir|bulk
|+|Fe bcc Green's function Radial Noco spin offdiagonal|greens/GreensFunction_mperp_ydir|bulk
|+|GdCu Green's function interorbital elements|greens/GreensFunction_InterOrbital|bulk
|+|Greens Function intersite single shell|greens/GreensFunction_IntersiteSingleShell|bulk
||Greens Function intersite single shell|greens/GreensFunction_IntersiteGammaNoGamma|bulk
|+|Greens Function intersite multiple shells|greens/GreensFunction_IntersiteMultipleShells|bulk
|+|Greens Function intersite shell construction|greens/GreensFunction_IntersiteShellConstruction|bulk
|+|Greens Function intersite shell construction|greens/GreensFunction_IntersiteShellConstructionFilm|bulk
