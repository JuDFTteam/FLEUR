# WhenDoINeedSOC

# When do I need SOC? 

*   (a) When calculating heavy atoms ( 6th or 7th period, sometimes also 5th) relativistic effects are important. Therefore, here it is often important to include SOC. Fine details of the bandstructure (e.g. heavy-hole / light-hole splitting in semiconductors) of course are also affected in lighter elements. 
*   (b) Calculating the magneto-crystalline anisotropy (MCA) for all magnetic materials requires the inclusion of SOC. Now, also the relative orientation of the axis of magnetization to the lattice has to be specified. 
*   (c) Whenever knowledge of the orbital moment is required. 
*   (d) When the magnetic structure is influenced by SOC effects (Dzyaloshinskii-Moriya interaction).

# HowDoICalculateWithSOC

# How do I calculate with SOC? 

This section describes the "usual" treatment of spin-orbit coupling for collinear calculations, where we use the so-called 2nd variation, i.e. the SOC operator is introduced in the subspace of eigenstates obtained without SOC (cf. C. Li et al.,[![Symbol - externer Link][2]Phys. Rev. B **42** 5433 (1990)][2]). For a treatment in 1st variation, non-collinear calculations can be performed. 

(a) Recompile with -DCPP_SOC. Make 'rminv' to be sure, all relevant .o-files are removed prior to compilation. 

(b) In the inp-file, set l_soc true in the line of the input file that is shown below. 



    |   0.00000   0.00000,l_soc=F,spav=F,off=F | 
    

Click here for [input][2] file info. 

In magnetic calculations, you have to specify two angles that give the relative orientation of the magnetisation axis to the lattice. The two numbers preceding 'l_soc' are the polar angle Î¸ and the azimuthal angle Ï† specifying this orientation. E.g. if the magnetic moments should point in z-direction, Î¸=0.0 and Ï†=0.0. 

The x-direction would be specified as: 



    |   1.57080   0.00000,l_soc=T|      ( &theta; = &pi;/2 )
    

and the y-direction: 



    |   1.57080   1.57080,l_soc=T|      (&theta; and &phi; = \pi/2 )
    

At each angle, you can replace the last two digits with 'pi' or 'dg' to indicate that you enter the number in multiples of Ï€ or in degrees, thus for `(&theta; and &phi; = &pi;/2 )` you can also write 

    |   0.5  pi   90.  dg,l_soc=T|
    

(c) Increase the number of states that are calculated per k-point. SOC is included in 2nd variation (see above) so enough unoccupied states have to be calculated in 'first variation' (i.e. the unperturbed problem = without SOC). In the parameter file ([fl7para][3]) increase neigd: 



    |Number of (occupied) bands    |
     |parameter (nobd=115,neigd=115)|
    

Then we have to adjust the size of the energy window: 



    |Window # 1                    |
     |  -0.50000   0.80000 140.00000|
                   ^^^^^^^
                   this should be set to, e.g. 1.8 htr.
    

Generally, it is recommended to check the convergence of the result with the number of additional states. 



## Special features:

You can use further switches in lines 24 and 31 of the inp file: 

    24|lpr=0,form66=F,l_f=F,eonly=F,eig66=T,soc66=T  |
     31|   0.00000   0.00000,l_soc=T,spav=T,off=T,1011|
    



*   If you set `spav=T`, the spin-orbit operator is constructed from a spin-averaged potential. 
*   If you set `off=T`, only the spin-orbit contribution of certain muffin tins is considered. These muffin tins are specified with a binary number at the end of line 31 (in the present example, spin-orbit coupling is considered for the 1st,3rd,4th atom type but ignored for the 2nd atom type). 
*   If you set `eig66=T`, then `soc66` determines whether the 'eig' file [that is created and read in different runs of the program] contains the eigenvectors from the first variation (without spin-orbit coupling, `soc66=T`) or from the second variation (`soc66=F`).

 []: http://dx.doi.org/10.1103/PhysRevB.42.5433
 
# SOC-RestrictionsAndProblems

Note, that using SOC is incompatible with the following features of the program: 



1.  Eigenvalue-parallelization. 
2.  Second-variation, Wu-diagonalization. 
3.  Spin spiral calculations 
4.  Constrained-moment calculations 
5.  Relaxation of mag. moments 

For (1) it is recommended to use SOC in 1st variation (as part of a non collinear calculation) 

In second variation with k-point parallelization use a number of k-points that is a multiple of the number of processors you use. 

For (3) you can add SOC in a [perturbative approach][1]. 

Furthermore, we note that for using LO's and SOC the COMPLEX version of the code is recommended, especially when forces are required. 



### Common problems: 

*   Calculation of the magneto-crystalline anisotropy energy (MAE) requires good convergence with respect to the number of k-points, the number of additional states, the 'extra vacuum' d-tilde in film-calculation etc. Normally, MAE is a very tiny quantity so it is very sensitive to all kinds of cut-offs. 

*   In calculating DOS or partial charges, don't be confused by the fact that SOC doubles the number of states in each spin-channel: this is a projection of the spins on the 'usual' spin - quantization axis. Note, that this also doubles the "sum of valence eigenvalues" sometimes used for MAE calculation with the force theorem.

 
