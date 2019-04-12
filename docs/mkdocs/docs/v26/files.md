FLEUR input- & output-file
======================
(Please note, that not all documentation covers v0.27)

**This you have to specify** 

* [[input file for the input generator]]
or
* The [inp](../inpfile.md) file (deprecated in v0.27)


## These you have to specify sometimes  

* The [nocoinp](#nocoinp) file
* The [qpts](#qpts) file
* The [plotin-file](#plotin) (deprecated in v0.27)
* The [band_inp](#band_inp) file
* Other [small useful](#small-usefull-files) files
* [Wannier-related input and output files](wannier.md)

## These you can modify 

* The [enpara](#enpara) file  (deprecated in v0.27)
* The [kpts](#kpts) file (deprecated in v0.27, replaced by list in inp.xml)
* The [sym.out](#sym.out) file
* The [fl7para](#fl7para) file (deleted in v0.27)
* The [plot_inp](#plotin) file

## Output produced by FLEUR 

* The [out](#out) file
* The [DOS.x](#dosx) files
* The [dosinp](#dosinp) file (deprecated in v0.27)
* The [Charge_and_Potential](#charge-and-potential-files) files
* The files produced in a [Jij calculation](#jij)
* Other [strange_unformatted](#strange-files) files
* Additional files produced in a [Hybrid functional](#hybrid) calculation

#nocoinp
The nocoinp file contains the settings needed for a calculation using Non-Collinear Magnetism

Note: This is a fixed-format file, so adhere to the format given below!

```
atom-type   1,l_relax=F
alpha =  0.0000000000,b_cons_x =  0.0000000000
beta  =  1.5707963268,b_cons_y =  0.0000000000
                                  0.0000000000
**********logical parameters******************
l_ss=F,l_mperp=F,l_constr=F,l_disp=F,sso_opt=FFT
mix_b= 0.500
qss=(  1.0000000000,  0.5000000000,  0.0000000000)
qsc=( 10.0000000000,  1.0000000000,  1.0000000000)
```

1. Select the atom type	atom-type: clear. l_relax: do you want to relax the direction of the moment localised at atom 1?
2. Angles	alpha: 1st angle that determines magnetic strucuture. Equal to "phi" in spherical coordinates. b_cons_x: Constraint-Field in x-direction. Is determined self-consistently if l_constr=T.
3. Angles beta: "theta" (measured from the z axis) in spherical coordinates.
For alpha and beta (lines 2-3) you can replace the last two digits by 'pi' or 'dg' to indicate that you enter the angles in multiples of Pi or in degrees, thus 'beta = 1.5707963268' is equivalent to 'beta = 0.50000000pi' or 'beta = 90.0000000dg'.
4. Reserved for future use.
For each atom-type these four lines have to be supplied!

5. empty line. Now the atomtype-dependent information is finished.
6. Logical switches	l_ss=T/F : Spinspiral, if l_ss=T then enter the spiral-vector in line 8. l_mperp=T/F : Do you want output of the magnetisation perpendicular to chosen axis (determined by alpha and beta)? l_constr=T/F : Do we constrain the moments or not. So far l_relax and l_constr exclude each other. sso_opt=??? : Three logical switches for (spin spirals + spin orbit).
7. Mixing mix_b= : Mixing-factor; if l_constr=T then mixing of Constraint-field in this case mix_b= 0.5 should work fine in case of l_relax=T Mixing of input/output-direction of moments you can choose mix_b>1 (e.g. 4)
8. Spin spiral vector qss : Measured in reciprocal lattice vectors. 
To each atom with a basis vector τ this will add an angle 2 π (q . τ) to α defined above (line 2).
9. [optional] denominator of spin spiral vector:
if this line (starting with 'qsc=') is provided, the components of qss are divided by the numbers given here
When you want to calculate Heisenberg interaction parameters Jij, you'll need additional switches. This is how nocoinp file looks in this case:

```
atom-type   1,l_relax=F,l_magn=T,M=2.96337
alpha =  0.0000000000,b_cons_x =  0.0000000000
beta  =  1.5707963268,b_cons_y =  0.0000000000
                                 0.0000000000
**********logical parameters******************
l_ss=T,l_mperp=F,l_constr=F,l_disp=F
mix_b= 0.500,thetaJ=  0.5235987756,nsh=  50
qss=(  0.1000000000,  0.5000000000,  0.0000000000)
```

1. Describe the atom-type l_magn: Is this atom-type magnetic?
         M: What is the value of its magnetic moment (don't forget the sign!)?
3. Warning:
   Angle beta will be taken into account only if l_disp=T (line 6),
   otherwise it is replaced with the angle thetaJ (line 7)
6. The dispersion switch l_disp: If l_disp=T, Force theorem is used to calculate the sum of
           eigenvalues for each qss predefined in qpts file. In this case,
           Jij parameters are not calculated.
7. Cone angle
  thetaJ: The angle used as beta (line 3) in Jij calculations
               nsh: How many shells of neighbours do you want to consider?
8. Warning: When l_J=T in the input file, qss defined here is ignored; qss vectors from qpts file are used instead.

# qpts
The qpts file contains a list of all spin-spiral vectors for a calculation of [[Heisenberg interaction parameters Jij]].
It has the following general format:

```                                                                                                                                         
  281
   .0000000000    .0000000000    .0000000000
   .4642857143    .4642857143    .4642857143
   .3928571429    .4642857143    .4642857143
 ...
```                                                                                                                                         
Where

* The first line specifies the total number of qss vectors
* A list of qss vectors follows, where their coordinates are measured in reciprocal lattice vectors

#plotin 
The plotin-file is the input file for the old charge-density 
and potential plotting subroutine. If you do not have no good reasons to stay with this
input format it is suggested that you use the new [plot_inp](#plot_inp) file instead.

Example of a file named "plotin":
```
 2-Dim=T
 Zahl   1
 xz_cut          test
  0.000000  0.000000 -19.000000
  1.000000  1.000000 -19.000000
  1.000000  1.000000  19.000000
 xPkt   100   yPkt   300
```
Creates a set of 100x300 charge density values
in a plane spanned by the vectors that connect the
3 positions (0,0,-19), (1,1,-19) and (1,1,19), 
where the first two coordinate values
are given in internal coordinates and the third
(z-value) is in atomic units.
The output file will be names zx_cut.

#band_inp

The band_inp file is a small input file used in the new bandstructure generation (see: [[Using the new FLEUR mode]]). It just contains a list of points in the Brillouine zone in the following format:

* the first character of the line is the Label of the point (use small letters for Greek symbols)
* the next three numbers are the coordinates of the point in internal coordinates. (For 2D band structures also three numbers have to be given but the last one will be ignored.) 

At the [Bilbao Crystallographic Server ](http://www.cryst.ehu.es/cryst/get_kvec.html) you find the ''k''-vectors types for all 230 space groups.

Example:
```
 g 0.0 0.0 0.0
 X 0.5 0.0 0.0
 M 0.5 0.5 0.0
 g 0.0 0.0 0.0
```

# Small useful files
In some cases, small files are used to influence the way,
how FLEUR  behaves during a run:

### qfix

Contains: a logical "T" or "F", to decide, how the charge is normalized, in case this has to be done. If "T" is set (or this file is not present), renormalize the charge everywhere in space, if "F" is read, renormalize only the interstitial charge.

Use: if the atomic positions have been changed, e.g. after a relaxation step, and the last converged charge density is used for the new calculation. 

### eps_force


Contains: a real, e.g. "   0.0005", to determine, when the convergence of the forces is good enough to make a new input file (inp_new). The default (if this file is not present) is 0.00001 htr/a.u.
Use: if the required precision is low, e.g. as long as you employ a steepest-descent algorithm for relaxation. 

### orbcomp

Contains: an integer, e.g. "3", to find out the atom, where an orbital decomposition should be made. If the DOS-flags are not set to `DOS=T` and `ndir=-3` this file is ignored, if the flags are set and no file is present, a layer-decomposition of the DOS is made.

### orbcomprot

**Contains: three lines, specifying the  Euler-angles {$\alpha, \beta$} and {$\gamma$} to rotate the coordinate frame in which the orbital decomposition (see above) is performed.

### apwefl

Contains: input parameters for applied external fields. It consists of one or two lines in free format, the first line contains a real number, specifying the location (in atomic units from the vacuum boundary) of the plane with the external charges as explained in the section on electric fields. The second (optional) line is used to specify the amount of charge that should be placed on this planes, one real number per vacuum (i.e. in case of inversion or z-reflection symmetry one number, otherwise two. For the syntax for inhomogeneous fields, read the page about electric fields.

### mfee

If present, an additional, external magnetic field is added. The syntax of the file is: One line per atom type; in each line the fixed-format "(i2,1x,f8.5)" is used, i.e. 2 digits for the atom type, space, and 8 digits for the value (in Hartree).

### n_mmp_rot

This file allows you to rotate the LDA+U density matrices by some angles (&theta;, &phi;) in real space. One line (free format) is reserved for each atom with an "U". The meaning of these angles is the same as in spin-orbit coupling mode, i.e. &theta; is a rotation around the y-axis, &phi; rotates afterwards around z. 

Use: e.g. to stabilize in an open d-shell an orbital moment in arbitrary directions, it is convenient to specify the density matrix with the moment along z (more or less diagonal) and then rotate in the required direction.  

### vca.in

Necessary for the "virtual crystal approximation", i.e. to tune the nuclear number between two integer values to simulate an alloy. E.g. a "element" with Z=25.5 could mean a Mn'_x_'Fe'_(1-x)_' alloy with x=0.5.  The file consists of lines with an integer (for the atom-type) and a float (added or subtracted charge at this atom). E.g. "1 0.5" means that half an nuclear charge should be added at the first atom type. The electronic charge should then be adjusted accordingly to ensure charge neutrality.

# enpara

Energy parameters for the linearized radial basisfunction 

There are, in principle, three ways to specify the energy parameters. These types are

* The normal "energy" format
* The "simple" n-format
* The floating energy parameters

### The normal "energy" format 

For each atom type in the input file the linearization energies are given for
all l-values.
Additionally you can specify:

* A mixing factor for a simple mixing scheme of these parameters 
* Logical switches allowing to fix some parameters
* The skiplo value specifying the number of states considered to be covered by LO's  

Depending on your setup you will also find:

* Parameters for both spins
* Energy parameters for the z-dependend basis in the vacuum for a film calculation
* Energies for local orbitals

Example (film, 1 atom type, 2 spins, no LO's):
```
     energy parameters for window 1 spin 1 mix=  1.000000
     atom     s        p        d        f
 --> 1   -0.31794 -0.25299 -0.21364 -0.21967 change: TTTT skiplo:   0
  vacuum parameter= -0.22066 change: T second vacuum= -0.22066
     energy parameters for window 1 spin 2 mix=  1.000000
     atom     s        p        d        f
 --> 1   -0.29891 -0.21597 -0.19182 -0.18453 change: TTTT skiplo:   0
  vacuum parameter= -0.20192 change: T second vacuum= -0.20192
```
Another example (film, 2 atom types, 1 spin, 2 LO's):
```
     energy parameters for window 1 spin 1 mix=  1.000000
     atom     s        p        d        f
 --> 1   -0.26719 -0.21434 -0.23628 -0.20953 change: TTTT skiplo:   4
 --> lo  -2.20113 -1.34013
 --> change   T        T
 --> 2   -0.25728 -0.22454 -0.21551 -0.22500 change: TTTT skiplo:   4
 --> lo  -2.21847 -1.35312
 --> change   T        T
  vacuum parameter= -0.22516 change: T second vacuum= -0.22516
```

Notice, that the LO's are one "s" and one "p" orbital, so that for the determination of the valence energy parameters (first line) in total the lowest four energy levels per atom have to be skipped ("skiplo: 4").

### The simple n-format 

Sometimes, it is tedious to find out the correct ranges for the energy parameters and how many levels have to be skipped for the LO's. Then it might be more convenient to simply give the principle quantum number of the state, you want to treat in your valence window and what should be added as local orbital. One example, how this can be done is shown here:
```
     energy parameters for window 1 spin 1 mix=  1.000000
     atom     s        p        d        f
 --> 1    4.00000  4.00000  3.00000  4.00000 change: FFFF skiplo:   4
 --> lo   3.00000  3.00000
 --> change   F        F
```

This would be the typical setup for an element on the left side in the 3d-row, with 4s, 4p, and 3d orbitals (plus 4f, also this has to be specified), and two local orbitals for 3s and 3p. For example Sc would require this kind of setup, maybe also Ti. The integer values of all energy parameters for an atom type indicates, that all parameters for this atom should be determined from an atomic calculation with the actual potential. Also the input-generator sets up this kind of files (the "change" flags are set to "F" to keep this format).

Which energy parameters are in the end actually taken for the calculation, you can read from the out file. For the above example (Ti) the program calculates e.g.:
```
 Atom  1 4s branch from -1.78 to  1.61 htr. ; e_l =  0.2179
 Atom  1 4p branch from -0.89 to  2.06 htr. ; e_l =  0.3310
 Atom  1 3d branch from -9.99 to  0.51 htr. ; e_l =  0.3247
 Atom  1 4f branch from -9.99 to  2.64 htr. ; e_l =  0.6168
 Atom  1 3s branch from-19.22 to -1.75 htr. ; e_l = -1.7786
 Atom  1 3p branch from-15.82 to -0.86 htr. ; e_l = -0.8986
```
this corresponds to an enpara-file in the "energy format":
```
     energy parameters for window 1 spin 1 mix=  1.000000
     atom     s        p        d        f
 --> 1    0.21790  0.33100  0.32470  0.61680 change: TTTT skiplo:   4
 --> lo  -1.77860 -0.89860
 --> change   T        T
```
The program will set the enpara-file back to this format, if you change the F's to T's, to allow an update of the enpara-file.

Please note, that this format causes problems if it is used with very light atoms, e.g. H, since the 4f-parameter cannot be determined (for some MT-radii, no bound state can be found).

### The floating energy parameters 

Another possibility to specify the energy parameters is, to reference them to a certain value of the radial potential in the muffin-tin (at 1/4 of the distance to the MT-radius). While this has of course the disadvantage, that energy parameters cannot be simply transfered from one   muffin-tin radius to the other, the energy parameters are more stable with respect to strong fluctuations of the potential, as occurs e.g. at the beginning of the SCF cycle. This can be toggled by setting a switch in the [[input]]-file (inp, line 26 in the example) from 0 to 1. An example (like above, for Ti) would be:
```
     energy parameters for window 1 spin 1 mix=  0.300000
     atom     s        p        d        f
 --> 1    9.10093  9.15589  9.17100  9.18091 change: TTTT skiplo:   4
 --> lo   7.13609  8.00391
 --> change   T        T
```
The reference value (potential in the muffin-tin) can be found in the out-file:
```
 Reference energies for energy parameters
 spin 1, atom type  1 =   -8.870484   r= 0.64167
```
adding this value to the  floating energy parameters gives you back the normal "energy format", very similar to the values of the last section.



#kpts

The kpts-file contains a list of all k-points. You have to specify a kpts-file or use the in-build k-point generator.

It has the following general format:
```
  2      140.0000000000
  13.00000  14.00000  15.00000   8.00000
  ....
```

Where

* the first line specifies the total number of k-points (here 2) and a scaling factor (here 140.0).
* the following lines give the (x,y,z) values of the k-point (here 13.0/140.0, 14.0/140.0 15.0/140.0) and its relative weight. The total weight is simply obtained by summing up all the relative weights.

### Notes:

* Instead of creating the file by yourself, FLEUR can generate it for you, using the parameters in the [[input]] file.
* If you switch to another k-point set, make sure to use the proper value for nkptd in file [[fl7para]], or delete it in order to make FLEUR recreate it with default values.
* In case of a 2D-calculation for a [[film-setup]]. Only the (x,y) coordinates are given and the weight is in place of the z-column
* If you set tria=T (which is very useful for obtaining good [[Density of states]]) the [[k-point generator]] will generate additional information at the end of the file.

# sym.out

Despite its name, this is an input file for FLEUR. It contains the symmetry information generated by the inp- file generator (hence its name :-) ). You probably should never edit it!

# fl7para

The file fl7para is read in the very beginning of a FLEUR run and contains the dimension parameters of the arrays that are to be allocated. Each value is paired with the variable name in FLEUR, and preceded by a line of brief explanation. If not existent at program startup, it will be created automatically with reasonable default values, a lot of which are constructed out of the statements in the input file.

Note that newer versions of Fleur do not automatically generate an fl7para file. If you want to modify values controlled by the  fl7para file, you can create a fl7para file using the template from the out file.

### Notes 

* It is a fixed format file, edit with care.
* In many cases you do not need to touch this file.
* One of the other cases is when you need more unoccupied states: increase neigd then.
* Other examples include the case when you increase values in the input file. Then you need to increase to proper dimensions in fl7para accordingly.

# plot_inp

The plot_inp file defines the parameters used to generate a charge density plot (parameter `iplot=t` in the inp file).

Example file:
```
 2,xsf=t
  &PLOT twodim=t,cartesian=t
   vec1(1)=10.0 vec2(2)=10.0
   filename='plot1' 
  /
  &PLOT twodim=f,cartesian=f
   vec1(1)=1.0 vec1(2)=0.0 vec1(3)=0.0
   vec2(1)=0.0 vec2(2)=1.0 vec2(3)=0.0
   vec3(1)=0.0 vec3(2)=0.0 vec3(3)=1.0
   grid(1)=30  grid(2)=30  grid(3)=30
   zero(1)=0.0 zero(2)=0.0 zero(3)=0.5
   filename ='plot2' 
  /
```

The first line specifies the no of plots to generate (2) and the output format. Setting `xsf` to true will generate output for the XCrysDen visualisation program.

The parameters of each plot are given in a namelist input below. 
The following set of input variables can be specified:

* `twodim=t/f`: twodimensional or threedimensional plot.
* `cartesian=t/f`: are the vectors given in cartesian or internal coordinates.
* `vec1`,`vec2`,`vev3`: The vectors spanning the plotting plane/volume. All components of these vectors default to 0.0, so that you have to give only the non-zero elements. For a 2D plot only `vec1` and `vec2` are used.
* `zero`: A vector shifting the origin of the plot volume, i.e. the corner from which vec1,vec2,vec3 are measured (defaults to (0.0,0.0,0.0)).
* `grid`: A integer vector specifying the grid size (two/three values for 2D/3D plots). All grid sizes default to 100.
* filename: A name for the file to store the data, or for the dataset in the xsf file.

There are additional options concerning '''noncollinear magnetization''':
```
 1,xsf=t,polar=t
  &PLOT phi0=0.5 unwind=t
  ...
```
* `polar=t/f`: In the case `polar=t`, additional files are created that allow to plot the magnetization in terms of absolute value ('mabs') and angles theta,phi ('mtha','mphi'). The default is `f`, the statement `polar=f` can be skipped. Note, that `polar` is not part of the namelist input.
* `phi0`: In the case `polar=t`, `phi0` specifies the range of the azimuth angle phi: phi is given in the interval [`phi0`*pi-pi, `phi0`*pi+pi]. The default is `phi0=`0.0.
* `unwind=t/f`: Relevant for spin-spiral plots. In the case `unwind=t`, the plot shows {$ {\rm m}_{\rm plot}({\rm r})= {\rm R}_{-{\rm q} . {\rm r}}\,{\rm m}({\rm r}) $}, i.e. the plot shows the magnetization rotated by the angle -q.r. In this case, the plotted {$ {\rm m}_{\rm plot} $} is periodic in the computational unit cell.

# out


The `out` file is the basic output file generated by the `fleur`-code. \
It contains a huge amount of information, normally it is not necessary to \
go through this file in detail, but in case of errors and problems it is in \
most cases useful to refer to this file.

The overall layout of this file is:

* a repetition of the input parameters
* symmetry information (stars & lattice harmonics)
* for each iteration:
  * multipoles, potential-density integrals
  * wavefunction parameters
  * eigenvalues for each k-point
  * Fermi level
  * partial charges & energy parameters
  * core-levels
  * magnetic moments, if any
  * input-output charge density distance 

If you compile `fleur.x` with `-DCPP_HTML`, you get a more structures, html-formated version of this file.

When you set up the starting density, you will also obtain atomic
levels of your constituent atoms in the out-file, e.g.
```
   n  kappa  l    j       occ.   eigenvalue (har)  <r>

 spin No. 1
   1   -1    0   0.5      2.00 -3324.468879    0.015793
   2   -1    0   0.5      2.00  -596.793964    0.065786
   2    1    1   0.5      2.00  -573.528807    0.053891
   2   -2    1   1.5      4.00  -488.314248    0.062935
   3   -1    0   0.5      2.00  -143.977036    0.171090
   3    1    1   0.5      2.00  -133.441852    0.161585
   3   -2    1   1.5      4.00  -114.369584    0.178084
   3    2    2   1.5      4.00   -97.056332    0.154532
   3   -3    2   2.5      6.00   -93.010528    0.159518
   4   -1    0   0.5      2.00   -32.992819    0.377147
   4    1    1   0.5      2.00   -28.414446    0.377905
   4   -2    1   1.5      4.00   -23.706774    0.412847
   4    2    2   1.5      4.00   -16.248825    0.416371
   4   -3    2   2.5      6.00   -15.370478    0.427488
   4    3    3   2.5      6.00    -5.679662    0.438154
   4   -4    3   3.5      8.00    -5.477905    0.444295
   5   -1    0   0.5      2.00    -5.745394    0.829854
   5    1    1   0.5      2.00    -4.149632    0.883427
   5   -2    1   1.5      4.00    -3.230129    0.974662
   5    2    2   1.5      4.00    -1.013890    1.215674
   5   -3    2   2.5      6.00    -0.903410    1.261554
   6   -1    0   0.5      2.00    -0.512439    2.182490
   6    1    1   0.5      2.00    -0.205150    2.715903
   6   -2    1   1.5      1.00    -0.132442    3.073088
```
(here for Bi). This might help to see, which states are well localized
(low energy and small radial extension <r>, in this case up to the 5p
levels), which are clearly valence states (the last three one: 6s, 6p1/2 and
6p3/2) and semicore states (5d3/2 and 5d5/2), which could be included
in the valence window using  local orbitals.


