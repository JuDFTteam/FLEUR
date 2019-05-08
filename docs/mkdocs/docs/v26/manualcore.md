# ManualCore-levelSetup

There are two cases when the use this feature might be necessary: 

*   when a unusual setup of the core-levels is wanted, e.g. for the 4f's in Lanthanides 
*   when there is no default core-level setup, as in the Actinides (see below) 



### Example Tb

If different setup of the core-levels than the default in the program is required, it is possible to specify a file named corelevels.Z, where Z is the nuclear number of the atom corresponding to the core-level setup given in the file. The format is as follows: 



    19
      1 -1  1  1  (1s)
      2 -1  1  1  (2s)
      2  1  1  1  (2p 1/2)
      2 -2  2  2  (2p 3/2)
      3 -1  1  1  (3s)
      3  1  1  1  (3p 1/2)
      3 -2  2  2  (3p 3/2)
      3  2  2  2  (3d 3/2)
      3 -3  3  3  (3d 5/2)
      4 -1  1  1  (4s)
      4  1  1  1  (4p 1/2)
      4 -2  2  2  (4p 3/2)
      4  2  2  2  (4d 3/2)
      4 -3  3  3  (4d 5/2)
      4  3  3  1  (4f 5/2)
      4 -4  4  0  (4f 7/2)
      5 -1  1  1  (5s)
      5  1  1  1  (5p 1/2)
      5 -2  2  2  (5p 3/2)
    

The first line gives the number of (relativistic) levels, the next lines specify 'n', 'kappa' and the occupancy of spin-up and -down of that level. In the above example, the 8 4f-levels of Tb have been put in the core. Note, that in this case the CPP_CORE compiler switch has to be set, otherwise the core-levels with energies above the interstitial average will leak out of the core potential. Moreover if a level is moved from the core to the valence region, or viceversa, the file 'enpara' which contains the energy of the local orbitals must be changed accordingly, as well as the value of the variables 'nlo' and 'llo' in the 'inp' file, which define the number and type of orbitals for which a local orbital is added. 



### Example Np

If you dare to work with very heavy elements, starting from Z=92, a careful setup of the core-levels is required. As an example take Np (Z=93), suppose you want to start with a configuration [Rn] 5f4, 6d1, 7s2. The corelevels.93-file looks like 



    27
      1 -1  1  1  (1s)
      2 -1  1  1  (2s)
      2  1  1  1  (2p 1/2)
      2 -2  2  2  (2p 3/2)
      3 -1  1  1  (3s)
      3  1  1  1  (3p 1/2)
      3 -2  2  2  (3p 3/2)
      3  2  2  2  (3d 3/2)
      3 -3  3  3  (3d 5/2)
      4 -1  1  1  (4s)
      4  1  1  1  (4p 1/2)
      4 -2  2  2  (4p 3/2)
      4  2  2  2  (4d 3/2)
      4 -3  3  3  (4d 5/2)
      4  3  3  3  (4f 5/2)
      4 -4  4  4  (4f 7/2)
      5 -1  1  1  (5s)
      5  1  1  1  (5p 1/2)
      5 -2  2  2  (5p 3/2)
      5  2  2  2  (5d 3/2)
      5 -3  3  3  (5d 5/2)
      6 -1  1  1  (6s)
      6  1  1  1  (6p 1/2)
      6 -2  2  2  (6p 3/2)
      7 -1  1  1  (7s)
      6  2  1  0  (6d 3/2)
      5  3  2  2  (5f 5/2)
    

As valence electrons we include 7s, 7p, 6d, and 5f, the 6p are included as local orbitals. Then, the enpara-file would look like 



    energy parameters for window 1 spin 1 mix=  1.000000
         atom     s        p        d        f
     --> 1    7.00000  7.00000  6.00000  5.00000 change: FFFF skiplo:   0
     --> lo   6.00000
     --> change   F
    

The corresponding inp-file (bcc Np for simplicity) starts like 



    strho=F,film=F,dos=F,isec1=99,ndir= 0,secvar=F
     Np bcc                                                                          
     any any ,invs=T,zrfs=F,invs2=F,jspins=1,l_noco=F,l_J=F
         -3.33130000     3.33130000     3.33130000
          3.33130000    -3.33130000     3.33130000
          3.33130000     3.33130000    -3.33130000    -3.33130000     1.00000000
     pz     non-relativi
    
       1
     **********************************
     Np  93   22   10  677  2.810000  0.020000
    
      1,force =F,nlo= 1,llo=  1
       0.000000  0.000000  0.000000  2.000000
     **********************************
      10.600000  8.800000
     vchk=F,cdinf=F,pot8=F,gw=0,numbands=  0
     lpr=0,form66=F,l_f=F,eonly=F,eig66=F,soc66=T
       8
       1  0
     Window # 1
       -0.80000   1.00000  13.00000
    

i.e. we include 22 core-levels (up to the 6s), and 13 valence electrons (6 for the p local orbitals, plus 4+1+2=7 lapw's).
