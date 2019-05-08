# Density Of States

To calculate the density of States different switches and modes are present in the inp-file: 



*   Most important, `dos` hast to bet set to .true. 

    strho=F,film=F,dos=T
    



*   The energy window, where the DOS is calculates is specified at the end of the inp-file: 

    emin_dos=  -0.50000,emax_dos=   0.50000,sig_dos=   0.00050
    

Note, that this is in Hartree units and not referred to some Fermi-level. The DOS is broadened by a Gaussian with the width ` sig_dos ` (also in htr.). 



*   A little bit nicer DOS' can be produced with the tetrahedron (3D) or triangular (2D) method. For films this is used whenever the k-point set allows it (corner points are included). For bulk, you have to create a k-point set with ` tria=T ` in the inp-file. 

*   Finally, there is the switch ` ndir ` in the first line of the inp-file. it is used trigger the following: 
    
    *   Producing raw-data for a DOS in the a dosinp-file
    *   Integrated DOS programm gives DOS.x files 
    *   Orbital decomposition of the DOS 

Note that you have to delete the fl7para and kpts files if you want to calculate the DOS.

# Producing Raw-data For A DOS

1.  Setting dos=t and ndir=0 in the inp-file switches into the DOS-mode of FLEUR. In this mode a dosinp-file is created which contains the eigenvalues and thier relative weights for all kpts. 
2.  Similarly, one obtains a vacDOS file containing the LDOS in the vacuum by setting vacdos=t. 
3.  By setting `dos=t` and `ndir>0` one obtains a dosinp-file with additional symmetry-information which is used for band-structure plotting

# IntegratedDOSProgramm

## total & partial DOS 

By setting dos=t and ndir<0 the integrated DOS program is used. Three additional parameters in the inp file have to be specified in this case. The two energies emin\_dos and emax\_dos specify the energy-range in which the DOS is evaluated. The parameter sig_dos specifies the broadening, as shown below (where the last few lines of the input file are given. 



    vacdos=T,layers= 1,integ=F,star=T,nstars= 2
     50
     iplot=F,score=F,plpot=F
     0  0.000000  0.000000,nnne=  0,pallst=F
     xa=   2.00000,thetad= 300.00000,epsdisp=   0.00010,epsforce=   0.00010
     relax 000
     emin_dos=    -1.000,emax_dos=   0.00000,sig_dos=   0.00050
    

Three different negative values for ndir are implemented at the moment: 

*   ndir=-1: In this mode FLEUR calculates the DOS-output and stops without creating a new charge density. This is analogous to the old DOS-mode 
*   ndir=-2: In this mode FLEUR calculates the DOS-output and creates a new charge density afterwards. This is done only if it=itmax, so that one can use this setting to create a DOS during self-consistency iterations. 
*   ndir=-3: In this mode the Orbital decomposition of the DOS is calculated 

For information about the output see the documentation of the DOS.x file. 



## vacuum DOS 

If one wishes to calculate the DOS in selected parts of the vacuum region, the parameter ` vacdos=T ` has to be set. In combination with `dos=T` and `ndir=-1`, this will generally produce a file VACDOS.[1,2] which has a format as describec in the documentation of the DOS.x file. 

This is done by integrating over the 2-D unit cell and in z-direction one can either integrate (`integ=T`) up several layers (defined by `layers=` in the inp-file), or choose planes, where the DOS is evaluated. 



*   If `integ=F`, the line below should provide as many integer numbers, as the number of `layers` that has been defined. The z-position of a layer with number N is given by: 



    dvac    N
       z = ---- + --  (a.u.)
             2    10
    
    



i.e. in the sample inp the layer # 50 is at 11 a.u. or about 3.7 Angstroem above the topmost Cu atom. The format is (in FORTRAN) ` (20(i3,1x))`, i.e. a maximum of 20 planes can be defined. By default, the largest number (N) that can be entered is 250. 

Notice, that there is a parameter ` layerd ` in the fl7para file. Typically, this should be increased to the number of layers. 



*   If integration is choosen, then the next line in the inp-file defines pairs of z\_low and z\_up for integration. Again the pairsare given in integers as described above. The format here is ` (10(2(i3,1x),1x)) `, i.e. up to 10 regions can be integrated. 

*   If you are interested in the lateral modulation of the local DOS in the vacuum, you can plot the (symmetrized) Fourier-components of the DOS by setting `star=T`, as long as you are within 10 a.u. from the vacuum matching plane (i.e. dvac/2). The number of Fourier-components you want so see in your VACDOS.[1,2] file is defined by `nstars=`. 

*   Many more features are built in, but usually the ones described above should suffice. In any case they are mainly tested within collinear calculations (without SOC), so be careful.

# Orbital Decomposition Of The DOS

In cases it might be useful to obtain also the orbital decomposed density of states of a certain atom. To get a DOS.x file with this information, you have to specify 



    dos=T
    

and 



    ndir=-3
    

in the inp file. Additionally, a file named `orbcomp` has to be present in the working directory, where the number of an atom (counted in the order as it appears in the inp file) is specified. Then, the DOS-file contains information about the selected atom (not atom-type!) in the following form (the energy is in eV): 



    energy  total-DOS                     selected atom-type
                       s   p   p   p   d    d    d    d      d   f   f   f   f   ...
                           x   y   z   xy   yz   zx   x²-y²  z²  x³  y³  z³  x²y
    

(Note: you have to provide k-points in the whole Brillouin-zone to get a correct DOS here! Proper symmetrization will be implemented in future releases.) The other f-components can be found in the file ek_orco.[1,2], where also the orbital composition is given for each k-point. To get this file, 



    cdinf=T
    

should be set in the inp-file. 



### rotated coordinate frames

Since the decomposition in spherical harmonics is per default performed in the global coordinate frame, it can be sometimes useful to have a tool to perform the decomposition in a rotated coordinate frame. For this, another file `orbcomprot` has to be specified that contains three lines, specifying the Euler-angles \alpha, \beta and \gamma for the rotation. 

## layer resolved DOS

If a zero is specified in the `orbcomp` file, or no such file has been provided, then the DOS.x files will contain the information about the layer resolved DOS. This is especially meaningful for a film-calculation, where the film can be divided into layers; the boundaries of the layers are defined as (x,y)-planes that intersects halfway between two neighbouring atoms of different height in the film. The results are listed in the following form: 



    energy  total-DOS   layer            interstitial         muffin-tin 
                                    contribution to layer
                         1, 2, ... n      1, 2, ... n          1, 2, ... n
