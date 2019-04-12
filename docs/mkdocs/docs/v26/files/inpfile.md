# Input

In the following a short sketch on the old inp file is provided:

# Description of an example "inp" file 

Note that we use atomic units, i.e. the length is in Bohr-radii (1 a.u. = 0.5291772108(18) Å) and the energy is in Hartree: 1 htr = 2 Ry = 27.2113845(23) eV. For bulk or film calculations, the energy zero is the average interstitial potential or the vacuum energy, respectively. 

## The inp file 

      
    [01][11]|strho=T,film=T,dos=F,isec1=99,ndir= 0,secvar=F  
    [02][12]|Cu 3l    tests  
    [03][13]|squ p4m ,invs=T,zrfs=T,invs2=T,jspins=1,l_noco=F,l_J=F  
    [04][14]|   4.82381  
    [05][15]| 12.000000 15.000000  1.000000  
    [06][16]|rpbe   non-tivistic  
    [07][17]|igrd=1,lwb=F,ndvgrd=6,idsprs=0,chng= -.100D-11  
    [08][18]|iggachk=0,idsprs0=0,idsprsl=0,idsprsi=0,idsprsv=0  
    [09][19]| 2  
    [10][20]|**********************************  
    [11][21]|Cu  29    7    8  421  2.150000   .023  
    [12][22]|  
    [13][23]| 1,force =F,nlo= 0,llo=  
    [14][24]|   .000000   .000000   .000000  1.000000  
    [15][25]|**********************************  
    [16][26]|Cu  29    7    8  421  2.150000   .023  
    [17][27]|  
    [18][28]| 2,force =T,nlo= 0,llo=  
    [19][29]|  .500000   .500000  3.050000  1.000000  
    [20][30]|   .500000   .500000 -3.050000  1.000000  
    [21][31]|**********************************  
    [22][32]| 10.500000 10.000000  
    [23][33]|vchk=F,cdinf=F,pot8=F,gw=0,numbands=  0  
    [24][34]|lpr=0,form66=F,l_f=F,eonly=F,eig66=F,soc66=F  
    [25][35]|  8  8  
    [26][36]|  1  0  
    [27][37]|Window # 1  
    [28][38]|  -1.00000   0.20000  33.00000  
    [29][39]|   3.80000  
    [30][40]|gauss=F    .0020 tria=F  
    [31][41]|   0.00000   0.00000,l_soc=F,spav=F,off=F,01  
    [32][42]|frcor=F,slice=F,ctail=F,disp=F,kcrel=0,u2f=F,f2u=F,bmt=F  
    [33][43]|itmax= 8,maxiter= 19,imix= 7,alpha=  0.10,spinf=  1.00  
    [34][44]|swsp=F  0.00  0.00  
    [35][45]|lflip=F  1  1  
    [36][46]|vacdos=F,layers= 1,integ=F,star=F,nstars= 0     0.00     0.00     0.00     0.00,nstm=0,tworkf=  0.000000  
    [37][47]|  
    [38][48]|iplot=F,score=F,plpot=F,band=F  
    [39][49]|  0   .000000   .000000,nnne=  0,pallst=F  
    [40][50]|xa=   2.00000,thetad= 300.00000,epsdisp=    .00010,epsforce=    .00010  
    [41][51]|relax 000 001  
    [42][52]|emin_dos=  -0.50000,emax_dos=   0.50000,sig_dos=   0.01500  
    [43][53]|nkpt=   20  
    [44][54]|nqpt=  200  
    ALTERNATIVELY  
    [43][53]|nkpt=   20,nx=06,ny=06,nz=06,gamma=F  
    [44][54]|nqpt=  200,qx=06,qy=06,qz=06  
    

[(1)][55] 

    strho =[T,F] if true, a starting-density is generated
     film  =[T,F] selects film (T) or bulk (F) calculations
     dos   =[T,F] generate dos-output file dosinp and stops
     isec1 =[0-99] iterative diagonalization used after iteration# isec1
     ndir  =[0-5]  if dos=T and ndir>0, calculate symmetry information for 
                   bandstructures; indicates which symmetry operations to use.
                   In version 22o and higher this has new behaviour!!!!!!
     secvar=[T,F]  non-spherical Hamiltonian treated in second variation
    

[(2)][57] Title and/or comment line 

[(3)][58] 

    latnam=[squ,p-r,c-r,hex,hx3,...] selects type of lattice 
     spgrp =[p4m ,pmm ,cmm ,p3m1,...] selects space-group
     invs  =[T,F] T, if the system has inversion-symmetry  
     zrfs  =[T,F] T, if the system has z-reflection symmetry 
     invs2 =[T,F] T, if the vacuum planes have 2-dimensional inversion-s.
     jspins=[1,2] numer of spins: paramagnetic (1) or magnetic (2) calculation 
     l_noco=[T,F] T, if non-collinear calculation (prepare file 'nocoinp')
     l_J=[T,F] T for a calculation of Heisenberg Jij parameters (goes with l_noco=T) 
    

[(4)][59] in-plane lattice constant(s) (alternatively the bravais lattice matrix can be given) 

    a1 (,a2)
    

[(5)][60] c-axis 

    - in case of bulk -
      c-axis   lattice constant
      c-axis   lattice constant
      scale    scaling factor for all lattice constants & z-coordinates 
     - in case of film -
      dvac     vacuum boundary (for film=T, otherwise dvac=dtild)
      dtild    z-boundary for 3D-planewave box ( > dvac !)
      scale    scaling factor for all lattice constants & z-coordinates 
    

[(6)][62] Exchange Correlation Potential settings 

    xc-potential=[x-a, mjw, [![Symbol - externer Link][64]pz][64], [![Symbol - externer Link][65]bh][65], wign, hl, [![Symbol - externer Link][66]vwn][66], xal, [![Symbol - externer Link][67]l91][67], [![Symbol - externer Link][68]pw91][68], [![Symbol - externer Link][69]pbe][69], [![Symbol - externer Link][70]rpbe][70], [![Symbol - externer Link][71]Rpbe][71]]
     relativistic ... uses relativistic corrections of [![Symbol - externer Link][72]MacDonnald-Vosko][72]. 
    

Please note that relativistic corrections in conjunction with the GGA are currently not implemented. 

[(7)][72] 

    igrd  =[0,1]   igrd=0: no gradient correction 
     lwb   =[T,F]   use White & Bird trick (disabled)
     ndvgrd=[2,4,6] grid partition for calculation of derivatives
     idsprs=[0,1]   general GGA print-switch
     chng  =        lowest allowed density value to before stop
    

[(8)][73] various GGA print-switches 

[(9)][74] 

    ntype =[0-99] number of atom types
    

[(10)][75] separator 

[(11)][76] Atom 

    Name of atom type     [Va, H,...,Lw]
     Nuclear Number        [ 0, 1,...103]
     number of core levels [typically 1,3,7,...]
     l-expansion cutoff    [typically 6-12]
     muffin-tin gridpoints [odd number, typically > 301]
     muffin-tin radius     [choose non-overlapping]
     logarithmic increment [for normal radii & meshes 0.02 -0.03]
    

[(12)][77] input line for LDA+U 

[(13)][79] 

    number of equivalent atoms in this atom type
     force =[T,F]          calculate forces on this atom-type
     nlo   =[0-99]         number of local orbitals to use
     llo   =[0-99],...     l-values for local orbitals
    

[(14)][80] Positions 

    x,y,z coordinates of atom (x&y always in internal (relative) units, if film=F also z)
     scale    scales coordinates by 1/scale. If film=T, scales only x&y coordinates, if film=F also z
    

[(15-20)][81] same as [(10-14)][20] for the second atom type 

[(21)][82] separator 

[(22)][83] Planewave cutoff 

    gmax         cutoff for PW-expansion of potential & density  ( > 2*kmax)
     gmaxxc       cutoff for PW-expansion of XC-potential ( > 2*kmax, < gmax)
    

[(23)][84] logical switches 

    vchk    =[T,F]   check continuity of potential at muffin-tin & vacuum boundary
     cdinf   =[T,F]   calculates partial charges and continuity of density
     pot8    =[T,F]   if T, use potential from files pottot and potcoul
     gw      =[0,1,2] controls the ouptut for the GW code Spex 
     numbands= N      sets the maximal number of bands to N
    

[(24)][85] logical switches 

    lpr   =[0,1] if lpr.gt.0, then also list eigenvectors on output file
     form66=[T,F] gives a formatted eigenvector file (eig)
     l_f   =[T,F] calculate pulay-forces on atoms, otherwise only HF-force
                  master switch for Geometry optimizer
     eonly =[T,F] if T, no eigenvectors are dumped on file 'eig'
     eig66 =[T,F] if T: if 'eig' file exists use eigenvalues and -vectors from 'eig',
                        if 'eig' file does not exist create it and stop
     soc66 =[T,F] relevant only for spin-orbit calculations with eig66=T
    

[(25)][87] l-cutoffs for the non-spherical Hamiltonian for all atom-types. 

    Notice that this is assumed to be < 10.
    

[(26)][89] old switches (do use with maximum care) 

    number of windows [1,2 or more]
     lepr=[0,1] energyparameters given on absolute (0) or floating (1) scale
    

[(27)][90] separator for each window (lines 27 to 29 are repeated for each window) 

[(28)][91] Energy window 

    lower energy boundary for eigenvalues (in hartree units)
     upper  (these boundaries are not used on the Cray or if invs=F)
     number of electrons in the window
    

[(29)][92] cutoff for Plane wave expansion of wavefunctions 

    kmax determines basis size
    

[(30)][93] Fermi energy, k-integration, weights 

    gauss=[T,F] use gaussian smearing for calculation of fermi-energy & weights
                 if gauss=F & tria =F fermi smearing is used (recommended for self-consistency) 
     tkb  =      temparature for smearing with gauss or fermi-smearing method
     tria =[T,F] use triangular method (2D version of tetrahedron method).
    

[(31)][94] SOC switches 

    theta,phi   angles to specify the  spin-quantization axis if l_soc=T
     l_soc=[T,F] use spin-orbit coupling
     spav =[T,F] construct spin-orbit operator from spin-averaged potential
     off  =[T,F] only soc contributions from certain muffin tins are considered
                 (atom types are specified by binary number)
    

[(32)][95] more switches 

    frcor=[T,F] if T, use frozen core approximation
     slice=[T,F] if T, calculate a slice (parameters in line (39)
     ctail=[T,F] if T, make core-tail correction (reexpansion of core-tails)
     disp =[T,F] if T, calculate the distance of in- and output potential
     kcrel=[0,1] for 0 (1), a fully-relativistic (spin-polarized) core routine is used
     u2f  =[T,F] generates a formatted density/potential from unformated file f_unf
     f2u  =[T,F] generates unformatted files from formatted cdn_form
     bmt =[T,F] generates density 'cdnbmt' with magnetization in interst. and vac. set to zero
    

[(33)][98] Mixing 

    itmax  =[1-99]    number of iterations done in this run
     maxiter=[0-99]    number of iterations used for broyden-mixing 
     imix   =[0,3,5,7] type of mixing (straight, Broyden 1st and 2nd or Anderson)
     alpha  =[0.-0.99] mixing factor (if > 10.0, only mixing is performed)
     spinf  =[1-100.0] spin mixing factor enhancement
    

[(34)][99] Initial magnetic moments 

    swsp=[T,F]    if T, generate spin-polarized density from unpolarized 
     bmu's         moments of atom-types generated if swsp=T
    

[(35)][100] Flip spins 

    lflip=[T,F]   if T, flip spin-direction for selected atoms
     nflip=[-1,1]  flip spin-direction for atom-types where nflip=-1
    

[(36)][101] Layered vacuum DOS 

    vacdos=[T,F]  if T, in case of dos=T also the dos in the vacuum region is calculated
     layers=[0-99] number of layers, in which the vacuum dos is integrated (see next line)
     integ =[T,F]  if T, vacuum dos is integrated also in  z-direction
     star  =[T,F]  if T, star coefficients are calculated at values of izlay for 0th (=q) to nstars-1
     nstars=[0-99] number of star functions to be used (0th star is given by value of q=charge integrated in 2D)
     locx, locy: four real numbers that can be used to calculate local DOS at a certain vertical position z 
         (or integrated in z) within a restricted area of the 2D unit cell, the corners of this area is given 
         by locx and locy they are defined in internal coordinates
     nstm  =[0-2] 0: s-Tip, 1: p_z-Tip, 2: d_z^2-Tip (following Chen's derivative rule)
     tworkf= Workfunction of Tip (in hartree units) is needed for d_z^2-Orbital)
    

[(37)][102] 

    if integ=T this line defines the z_low and z_up for integration (in internal units)
     otherwise the z_values of the planes are entered.
    

[(38)][103] Charge/ Potential plotting 

    iplot=[T,F]   calculate a charge density plot
     score=[T,F]   if T, excludes the core-charge from the plot
     plpot=[T,F]   allows to plot the potential from potential-files
     band =[T,F]   simplifies the creation of band structure plots
    

[(39)][106] Charge density slicing 

    number of k-point which is used for a [slice][96] (k=0 : all k-points taken)
     lower boundary for eigenvalues in the slice 
     upper -"-
     nnne  = number of eigenvalue used for the slice (nnne=0 : all eigenvalues between boundaries taken)
     pallst=[T,F] set true if one plots states which lie above the fermi level
    

[(40)][107] Geometry optimizer 

    xa      = mixing parameter for geometry optimizer (2. or 3. is a good choice)
     thetad  = debye temperature used for first geometry optimization step
     epsdisp = if all displacements are < epsdisp, the program stops
     epsforce= f all forces  are < epsforce, the program stops 
    

[(41)][108] Geometry optimizer 

    relax: for each atom-types a triple of 0's or 1's specifies if the (x,y,z)
     coordinates can be relaxed; i.e. 001 means that only relaxation in z-direction is 
     allowed.
    

[(42)][109] DOS output parameter 

    emin_dos= set lower boundary of energy window of the DOS plot
     emax_dos= set upper boundary of energy window of the DOS plot
               (both values only affect the plot, not the energy window of eigenvalues specified above)
     sig_dos = Gaussian smearing factor used in the plot (if tetrahedron method is not used)
    

[(43)][110] k-points mesh (to be generated if not already existent) 

    nkpt      =  number of IBZ k-points to be generated (only very rough estimate; only used if nx/ny/nz not specified.)
     nx,ny,nz  =  x/y/z mesh for equidistant full BZ mesh to be generated. (This is optional, see above.)
     gamma     = [T,F] is a optional keyword; if true a k-point set will be generated, which includes the Gamma point as the  
                       first k-point
    

[(44)][111] Spin-spirals (qss) mesh, generated if not already existent 

    nqpt     = number of IBZ qss to be generated (used only if qx/qy/qz not specified)
      qx,qy,qz = x/y/z mesh for equidistant full BZ qss mesh to be generated (This is optional, see above).

 [11]: #r01
 [12]: #r02
 [13]: #r03
 [14]: #r04
 [15]: #r05
 [16]: #r06
 [17]: #r07
 [18]: #r08
 [19]: #r09
 [20]: #r10
 [21]: #r11
 [22]: #r12
 [23]: #r13
 [24]: #r14
 [25]: #r15
 [26]: #r16
 [27]: #r17
 [28]: #r18
 [29]: #r19
 [30]: #r20
 [31]: #r21
 [32]: #r22
 [33]: #r23
 [34]: #r24
 [35]: #r25
 [36]: #r26
 [37]: #r27
 [38]: #r28
 [39]: #r29
 [40]: #r30
 [41]: #r31
 [42]: #r32
 [43]: #r33
 [44]: #r34
 [45]: #r35
 [46]: #r36
 [47]: #r37
 [48]: #r38
 [49]: #r39
 [50]: #r40
 [51]: #r41
 [52]: #r42
 [53]: #r43
 [54]: #r44
 [55]: #l01
 [57]: #l02
 [58]: #l03
 [59]: #l04
 [60]: #l05
 [62]: #l06
 []: http://dx.doi.org/10.1103/PhysRevB.23.5048
 []: http://www.iop.org/EJ/abstract/0022-3719/5/13/012/
 []: http://www.nrcresearchpress.com/doi/abs/10.1139/p80-159
 []: http://dx.doi.org/10.1103/PhysRevB.45.13244
 []: http://dx.doi.org/10.1103/PhysRevB.46.6671
 []: http://link.aps.org/abstract/PRL/v77/p3865
 []: http://link.aps.org/abstract/PRL/v80/p890
 []: http://link.aps.org/abstract/PRB/v59/p7413
 []: http://dx.doi.org/10.1088/0022-3719/12/15/007
 [72]: #l07
 [73]: #l08
 [74]: #l09
 [75]: #l10
 [76]: #l11
 [77]: #l12
 [79]: #l13
 [80]: #l14
 [81]: #l15
 [82]: #l21
 [83]: #l22
 [84]: #l23
 [85]: #l24
 [87]: #l25
 [89]: #l26
 [90]: #l27
 [91]: #l28
 [92]: #l29
 [93]: #l30
 [94]: #l31
 [95]: #l32
 [98]: #l33
 [99]: #l34
 [100]: #l35
 [101]: #l36
 [102]: #l37
 [103]: #l38
 [106]: #l39
 [107]: #l40
 [108]: #l41
 [109]: #l42
 [110]: #l43
 [111]: #l44