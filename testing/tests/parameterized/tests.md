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

|     | Description                                  | directory name           | marks                 | Remarks            | CmdLine         | MPI |
| --- | -------------------------------------------- | ------------------------ | --------------------- | ------------------ | --------------- | --- |
| +   | Bulk Cu one iteration,bandstructure          | basic/CuBand             |                       |                    |                 | 2   |
| +   | Bulk Cu one iteration,ELPA                   | basic/CuBulk             | mpi,elpa              | Test for ELPA      | -diag,elpa      | 2   |
| +   | Bulk Cu one iteration,ScaLapack              | basic/CuBulk             | mpi                   | Test for ELPA      | -diag,scalapack | 2   |
| +   | Bulk Cu one iteration,Chase                  | basic/CuBulk             | mpi,chase             | Test for ELPA      | -diag,chase     | 2   |
| +   | Bulk Cu one iteration                        | basic/CuBulk             | fast,bulk             |                    |                 | 2   |
| +   | Bulk Cu one iteration,DOS                    | basic/CuDOS              | fast,bulk,dos         |                    |                 | 2   |
| +   | Bulk Cu one iteration,DOS,Orbital decomp.    | basic/CuOrb              | fast,bulk,dos,orbcomp |                    |                 | 2   |
|     | Bulk Co, DOS,MCD                             | basic/CoMCD              | bulk,dos,mcd          | MCD disabled       |                 | 2   |
| +   | Bulk Co, bandstructure, unfolding            | basic/CoUnfold           | band,bulk             |                    |                 | 2   |
| +   | Bulk Fe, Kerker preconditioner               | basic/Fe_Kerker          | bulk                  |                    |                 | 2   |
| +   | Bulk Fe fcc with relativistic core solver    | basic/Fe_fcc_kcrel       |                       |                    |                 | 2   |
| +   | Si with LOs                                  | basic/SiLO               | bulk                  |                    |                 | 2   |
| +   | Bulk PTO                                     | basic/PTO                | bulk                  |                    |                 | 2   |
| +   | Bulk PTO, SOC                                | basic/PTO-SOC            | bulk,soc              |                    |                 | 2   |
| +   | Bulk Fe, Tetrahedon method                   | basic/Fe_Tetra_noSYM     | bulk                  |                    |                 | 2   |
|     | LDA+U with AMF double counting and magnetism | basic/NiOldaUAMF         | bulk,ldau             | LDA+U AMF disabled |                 | 2   |
| +   | LDA+U with LF double counting and magnetism  | basic/NiOldaUFLL         | bulk,ldau             |                    |                 | 2   |
| +   | Crystal field output                         | basic/CrystalFieldOutput | bulk                  |                    |                 | 2   |

Testset: Films
----------------
|     | Description                    | directory name       | marks | Remarks                             | CmdLine | MPI |
| --- | ------------------------------ | -------------------- | ----- | ----------------------------------- | ------- | --- |
| +   | Fe Monolayer SOC               | film/Fe_1l_SOC       | soc   |                                     |         | 2   |
| +   | Fe Monolayer                   | film/Fe_1l           |       |                                     |         | 2   |
|     | Fe Monolayer Triangular method | film/Fe_1l_Tria      |       | not reproducible on INTEL compilers |         | 2   |
| +   | Si Film, plotting              | film/SiFilmPlot      |       |                                     |         | 2   |
| +   | Si Film, plotting,slicing      | film/SiFilmSlicePlot |       |                                     |         | 2   |
| +   | Pt 3 layers,soc,inversion      | film/Pt-3            | soc   |                                     |         | 2   |


Testset: Forces
------------
|     | Description                                  | directory name          | marks     | Remarks | CmdLine | MPI |
| --- | -------------------------------------------- | ----------------------- | --------- | ------- | ------- | --- |
| +   | Bulk GaAs, Relaxation, LDA+U                 | forces/GaAsMultiUForce  | bulk,ldau |         |         | 2   |
| +   | Bulk VO2, Relaxation                         | forces/VO2_forces       | bulk      |         |         | 2   |
| +   | Bulk VO2, Relaxation, different force levels | forces/VO2_force_levels | bulk      |         |         | 2   |
| +   | Bulk H2O, Relaxtion using BFGS               | forces/H2ORelaxBFGS     | bulk      |         |         | 2   |

Testset: DFPT
-------------------
|     | Description   | directory name              | marks | Remarks      | CmdLine | MPI |
| --- | ----------------------- | ----------------- | ----- | ------------ | ------- | --- |
| +   | Cu bulk fcc, Gamma only | dfpt/CuBulkGamma  | libxc | partly ready |         |   1 |
| +   | SiC bulk fcc, Gamma only BEC| dfpt/SiCBulk-BEC  | libxc |          |         |   1 |
| +   | Cu bulk fcc, K-mesh     | dfpt/CuBulkKmesh  | libxc |              |         |   1 |
| +   | C bulk fcc              | dfpt/CBulk        | libxc | partly ready |         |   1 |
| +   | C bulk fcc BEC             | dfpt/CBulk-BEC        | libxc | partly ready |         |   1 |
| +   | Graphene Film           | dfpt/GrapheneFilm | libxc |              |         |   1 |
| +   | V bulk bcc, mpi         | dfpt/VBulkMPI     | libxc |              |         |   2 |


Testset: Noco
----------
|     | Description                                   | directory name          | marks                        | Remarks           | CmdLine | MPI |
| --- | --------------------------------------------- | ----------------------- | ---------------------------- | ----------------- | ------- | --- |
| +   | Fe bct,noco, SOCL                             | noco/Fe_bct_SOC         | bulk,soc                     |                   |         | 2   |
| +   | Fe bct,noco, LOs                              | noco/Fe_bct_LO          | bulk                         |                   |         | 2   |
| +   | Fe bct,noco                                   | noco/Fe_bct             | bulk                         |                   |         | 2   |
| +   | Fe bct,noco,libxc                             | noco/Fe_bct_LibXC       | bulk,libxc                   |                   |         | 2   |
| +   | Noco, one atom, mag. in x direction           | noco/1atx               | bulk                         |                   |         | 2   |
| +   | Noco, one atom, mag. in y direction           | noco/1aty               | bulk                         |                   |         | 2   |
| +   | Noco, one atom, mag. in z direction           | noco/1atz               | bulk                         |                   |         | 2   |
| +   | Noco, one atom, mag. in non-sym direction     | noco/1at                | bulk                         |                   |         | 2   |
| +   | Noco,SOC, one atom, mag. in x direction       | noco/1atSOCx            | bulk,soc                     |                   |         | 2   |
| +   | Noco,SOC, one atom, mag. in y direction       | noco/1atSOCy            | bulk,soc                     |                   |         | 2   |
| +   | Noco,SOC, one atom, mag. in z direction       | noco/1atSOCz            | bulk,soc                     |                   |         | 2   |
| +   | Noco,SOC, one atom, mag. in non-sym direction | noco/1atSOC             | bulk,soc                     |                   |         | 2   |
| +   | Noco, SOC, two eq. atoms, relLOs              | noco/relLO              | bulk,soc                     |                   |         | 2   |
| +   | FFNNoco, one atom, mag. in x direction        | noco/1atFFNx            | bulk,hdf                     |                   |         | 2   |
| +   | FFNNoco, one atom, mag. in y direction        | noco/1atFFNy            | bulk,hdf                     |                   |         | 2   |
| +   | Fe fcc spin-spiral                            | noco/Fe_fcc             | bulk                         |                   |         | 2   |
| +   | Iron LO's and SOC test in FFN                 | noco/FeFFNLOsSOC        | bulk,soc,hdf                 |                   |         | 2   |
| +   | Fe monolayer fcc (110): SS                    | noco/FePt_film_SSFT     | film,spinspiral,forcetheorem |                   |         | 2   |
| +   | Fe monolayer fcc (110): SS +LO                | noco/FePt_film_SSFT_LO  | film,spinspiral,forcetheorem |                   |         | 2   |
| +   | Noco, Mn Monolayer SS q=1,0,0                 | noco/MnFilmSS           | film,spinspiral              |                   |         | 2   |
| +   | Noco, Mn Monolayer mag. in X                  | noco/MnFilmX            | film                         |                   |         | 2   |
| +   | Noco, Mn Monolayer mag. in Y                  | noco/MnFilmY            | film                         |                   |         | 2   |
| +   | Fe bct,noco,non-collinear,coretails           | noco/Fe_bct_ctail       | bulk                         |                   |         | 2   |
| +   | Noco, Mn Monolayer mag. in X, coretails       | noco/MnFilm_ctail       | film                         |                   |         | 2   |
| +   | Noco, one atom in x, noco IR starting density | noco/1atx_sdNocoIR      | bulk                         |                   |         | 2   |
|     | Fe bcc, Flipcdn and noco in MT,x-dir          | noco/Fe_bcc_FlipcdnXLDA | bulk                         | produces warnings |         | 2   |
|     | Fe bcc, Flipcdn and noco in MT,y-dir          | noco/Fe_bcc_FlipcdnYLDA | bulk                         | produces warnings |         | 2   |
| +   | relaxation feature of FFN in the MT           | noco/RelaxMT            | bulk,hdf                     |                   |         | 2   |

Testset: Experimental
----------
|     | Description                                    | directory name          | marks         | Remarks                   | CmdLine | MPI |
| --- | ---------------------------------------------- | ----------------------- | ------------- | ------------------------- | ------- | --- |
| +   | Test of the orbital polarization correction    | extra/Fe_bcc_OPC        | bulk,soc      |                           |         | 2   | 
| +   | Sourcefree magnetism and magnetization scaling | extra/Fe_bcc_SF_LDA     | bulk,hdf      |                           |         | 2   |
| +   | Bulk Cu one iteration                          | extra/CuBulkLibXC       | libxc,bulk    |                           |         | 2   |
| +   | Bulk Al one iteration, LibXC                   | extra/Al_libxc_PBE      | bulk,libxc    |                           |         | 2   |
|     | Test of GW interface 1                         | extra/gw1Interface      | bulk          | inp.xml files too old     |         | 2   |
|     | Test of GW interface 2                         | extra/gw2Interface      | bulk          | inp.xml files too old     |         | 2   |
|     | Sm jDOS decomposition                          | extra/SmAtomjDOS        | bulk,dos      |                           |         | 2   |
| +   | C: simple test for the Wannier code            | extra/Cwann             | bulk,wannier  |                           |         | 2   |
|     | TiO2 EELS spectrum                             | extra/TiO2eels          | bulk,eels     | inp.xml too old           |         | 2   |
| +   | Hubbard1 using SOC                             | extra/Gd_Hubbard1       | bulk,edsolver |                           |         | 2   |
| +   | Hubbard1 without sym                           | extra/Gd_Hubbard1_noSYM | bulk,edsolver |                           |         | 2   |
|     | diamond for one k-point with scan              | extra/Diamond_SCAN      | bulk,libxc    | SCAN has to be refactored |         | 2   |
| +   | 3D vector plots of the magnetization           | extra/PlotOnlyMT        | bulk,plot,hdf |                           |         | 2   |
| +   | density and potential plots, vector plots      | extra/PlotDenandPot     | bulk,plot,hdf |                           |         | 2   |


Testset: Hybrid-Functionals
----------
|     | Description          | directory name        | marks    | Remarks                                    | CmdLine | MPI |
| --- | -------------------- | --------------------- | -------- | ------------------------------------------ | ------- | --- |
| +   | GaAs PBE0            | hybrid/GaAsHybridPBE0 | bulk,hdf | runs too long                              |         | 2   |
| +   | Fe PBE0              | hybrid/FeHybridPBE0   | bulk,hdf | runs too long                              |         | 2   |
| +   | KCl PBE0             | hybrid/KClHybridPBE0  | bulk,hdf | runs too long                              |         | 2   |
|     | Mn Noinversion, PBE0 | hybrid/MnHybridNoinv  | bulk,hdf | runs too long, redundant to GaAs, Fe tests |         | 2   |



Testset: Greenfunctions
----------
|     | Description                                          | directory name                                       | marks | Remarks         | CmdLine | MPI |
| --- | ---------------------------------------------------- | ---------------------------------------------------- | ----- | --------------- | ------- | --- |
|     | Fe bcc Green's function                              | greens/Fe_bcc_GreensFunction                         | bulk  |                 |         | 2   |
|     | Fe Monolayer Green's function                        | greens/Fe_1l_GreensFunction                          | film  |                 |         | 2   |
|     | Greens Function MultiContour                         | greens/GreensFunction_MultiContour                   | bulk  |                 |         | 2   |
|     | Fe bcc Green's function Radial                       | greens/GreensFunctionRadial                          | bulk  |                 |         | 2   |
|     | Fe bcc Green's function Radial with local orbitals   | greens/GreensFunctionRadial_LO                       | bulk  |                 |         | 2   |
|     | Ho atom Green's function                             | greens/GreensFunction_HoAtom_SQA_theta               | bulk  |                 |         | 2   |
|     | Ho atom Green's function                             | greens/GreensFunction_HoAtom_SQA_phi                 | bulk  |                 |         | 2   |
|     | Ho atom Green's function                             | greens/GreensFunction_rotated_SQA_noco               | bulk  |                 |         | 2   |
|     | Fe bcc Green's function Radial Noco spin offdiagonal | greens/GreensFunction_mperp_xdir                     | bulk  |                 |         | 2   |
|     | Fe bcc Green's function Radial Noco spin offdiagonal | greens/GreensFunction_mperp_ydir                     | bulk  |                 |         | 2   |
|     | GdCu Green's function interorbital elements          | greens/GreensFunction_InterOrbital                   | bulk  |                 |         | 2   |
|     | Greens Function intersite single shell               | greens/GreensFunction_IntersiteSingleShell           | bulk  |                 |         | 2   |
|     | Greens Function intersite single shell               | greens/GreensFunction_IntersiteGamma                 | bulk  |                 |         | 2   |
|     | Greens Function intersite single shell               | greens/GreensFunction_IntersiteNoGamma               | bulk  |                 |         | 2   |
|     | Greens Function intersite multiple shells            | greens/GreensFunction_IntersiteMultipleShells        | bulk  |                 |         | 2   |
|     | Greens Function intersite shell construction         | greens/GreensFunction_IntersiteShellConstruction     | bulk  | takes too long! |         | 2   |
|     | Greens Function intersite shell construction         | greens/GreensFunction_IntersiteShellConstructionFilm | bulk  |                 |         | 2   | 

Testset: Wannier
------

|     | Description                              | directory name   | marks         | Remarks | CmdLine | MPI |
| --- | ---------------------------------------- | ---------------- | ------------- | ------- | ------- | --- |
| +   | Pt no-SOC, wannierlib total spread       | wannier/WannPt   | wannierlib,bulk  |         |         | 1   |
| +   | Pt SOC, wannierlib total spread              | wannier/WannPtSOC | wannierlib,bulk,soc |         |         | 1   |
| +   | fcc Fe FM noco, wannierlib total spread      | wannier/WannFeFM | wannierlib,bulk  |         |         | 1   |
| +   | fcc Fe AFM noco, wannierlib total spread     | wannier/WannFeAFM | wannierlib,bulk  |         |         | 1   |
| +   | fcc Fe AFM noco+SOC, wannierlib total spread | wannier/WannFeAFMSOC | wannierlib,bulk,soc |         |         | 1   |
| +   | Pt SOC, wannierlib real-space operators O(R) | wannier/WannPtSOCOps | wannierlib,bulk,soc |         |         | 1   |
| +   | bcc Fe FM collinear+SOC, wannierlib O(R)     | wannier/WannFeBccSOC | wannierlib,bulk,soc |         |         | 1   |
| +   | fcc Fe AFM collinear+SOC, wannierlib O(R)    | wannier/WannFeAFMColSOC | wannierlib,bulk,soc |         |         | 1   |
| +   | bcc Fe FM collinear no-SOC, wannierlib O(R)  | wannier/WannFeBcc | wannierlib,bulk |         |         | 1   |
| +   | fcc Fe AFM noco+SOC, wannierlib O(R)         | wannier/WannFeAFMSOCOps | wannierlib,bulk,soc |         |         | 1   |
