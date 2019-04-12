# HybridFunctionals

Hybrid functionals are not implemented for all features of FLEUR. In particular, the non-collinear calculations, the 1D & 2D version, band structures and the calculation of forces and relaxation are not implemented. Hence, hybrid functionals are not in the main release of FLEUR, so you have to use the source code provided below for all calculations associated with hybrid functionals. 


You can also take a look at [the talk about hybrid functionals][77] to get more information.

 [77]: https://www.flapw.de/pm/uploads/User-Documentation/schlipf.pdf



Introduction to Hybrid Functionals in FLAPW
=================
In hybrid functionals a certain fraction of non-local exchange replaces the corresponding amount of local exchange

\( E_{\mathrm{xc}} = E_{\mathrm{xc}}^{\mathrm{L}} + a \left( E_{\mathrm{x}}^{\mathrm{NL}} - E_{\mathrm{x}}^{\mathrm{L}} \right)\).

In *Fleur*, we have implemented the [PBE0](http://dx.doi.org/10.1063/1.472933) and the [HSE0](http://dx.doi.org/10.1063/1.2404663|HSE06) functional. In PBE0, the local part is the PBE-GGA functional and the mixing parameter is \(a = 25\%\). In HSE06, the local part we use the same mixing parameter \(a\) and also the PBE-GGA functional, but only the short-range component is mixed. The short-range is defined by the separation of the Coulomb potential

\( \frac{1}{r} = \underbrace{\frac{\mathrm{erfc}(\omega r)}{r}}_{\mathrm{short-range}} + \underbrace{\frac{\mathrm{erf}(\omega r)}{r}}_{\mathrm{long-range}}\)

with the complementary error function. The optimum choice for the adjustable parameter \(\omega\) was found to be \(0.11\).

In a generalized Kohn-Sham scheme, we calculate the additional contribution \(V^{\mathrm{NL}}_{\mathbf{G}\mathbf{G}'}\) of the non-local part via the double-integral

\( V^{\mathrm{NL}}_{\mathbf{G}\mathbf{G}'}(\mathbf{k}) = -\sum_{m}^{\mathrm{occ.}}\sum_{\mathbf{q}}^{\mathrm{BZ}} \int \mathrm{d}^3r \int \mathrm{d}^3r' \chi^*_{\mathbf{k+G}}(\mathbf{r})\phi^*_{m\mathbf{q}}(\mathbf{r'})v^{\mathrm{NL}}(|\mathbf{r}-\mathbf{r'}|)\phi_{m\mathbf{q}}(\mathbf{r})\chi_{\mathbf{k+G'}}(\mathbf{r'})\) ,

where the non-local potential \(v^{\mathrm{NL}}\) is either the Coulomb potential (PBE0) or the attenuated Coulomb potential (HSE06) and \(\chi\) are FLAPW basis functions. We can represent the basis functions \(\chi\) in terms of the eigenfunctions \(\phi\)

\( \chi_{\mathbf{k+G}}(\mathbf{r}) = \sum_{n}^{nbands} a_{\mathbf{k+G}}^{(n)} \phi_{n\mathbf{k}}(\mathbf{r})\),

where the \(a\)'s are the expansion coefficients. '''Note that the cutoff \(nbands\) can be chosen in the inp file and is the most important additional convergence parameter for hybrid functional calculations'''.

Implementation
-------
(More Details: [here](http://dx.doi.org/10.1103/PhysRevB.81.195117|here) and [here](http://dx.doi.org/10.1103/PhysRevB.84.125142|here))

After representing the basis functions via the eigenfuntions, we have to calculate an expression

\( \displaystyle V^{\mathrm{NL}}_{nn'}(\mathbf{k}) = -\sum_{m}^{\mathrm{occ.}}\sum_{\mathbf{q}}^{\mathrm{BZ}} \int \mathrm{d}^3r \int \mathrm{d}^3r' \phi^*_{n\mathbf{k}}(\mathbf{r})\phi^*_{m\mathbf{q}}(\mathbf{r'})v^{\mathrm{NL}}(|\mathbf{r}-\mathbf{r'}|)\phi_{m\mathbf{q}}(\mathbf{r})\phi_{n'\mathbf{k}}(\mathbf{r'}) \).

We calculate this introducing a new basis for the wavefunction products. The new basis consist of two types of functions, which are muffin-tin functions

\( M_{I\mathbf{k}}(\mathbf{r}) = \frac{1}{\sqrt N} \sum_{\mathbf T} R_I(|\mathbf{r} - \mathbf{T} - \mathbf{\tau}_I|) Y_I(\widehat{\mathbf{r} - \mathbf{T} - \mathbf{\tau}_I}) \exp[\mathrm{i}\mathbf{k}(\mathbf T + \mathbf{\tau}_I)]\)

and interstitial functions

\( M_{I\mathbf{k}}(\mathbf{r}) = \frac{1}{\sqrt{N\Omega}} \exp[\mathrm{i}(\mathbf{k+G}_I)\mathbf{r}]\Theta(\mathbf{r} \in \mathrm{IR})\).

Here \(N\) is the number of atoms in the crystal. \(\mathbf T\) is any crystal vector, \(\mathbf \tau_I\) is the position of the atom which is associated to the basis function \(I\), \(Y_I\) is a spherical harmonic where its subscript \(I\) is an abbreviation for \(Y_{l_Im_I}\), \(\Omega\) is the volume of the unit cell, and \(\Theta\) the Heaviside function which is only non-zero in the interstitial region (IR). As the basis consists of a mixture of two different basis functions, we refer to this basis as the mixed basis.

In the construction of the mixed basis more convergence parameters occur. The most important are \(gcutm\) a cutoff for the largest planewave used in the mixed basis and the atom specific \(lcutm\) which determines the maximum angular moment used. Note that as the mixed basis describes products of two wavefunctions the cutoff have to be larger. For \(lcutm\) twice the highest occupied angular moment of the atom is a reasonable choice (e.g. \(d\)-states \(\Rightarrow lcutm=4\), \(f\)-states \(\Rightarrow lcutm=6\)).

'''Nested Self-Consistency''':
In hybrid functional calculations the most expensive part (by a large margin) is the calculation of the non-local potential. Thus, we showed that it is very efficient to separate the self-consistency in two loops. In the outer loop, we calculate the non-local potential as well as the local potential. In the inner loop, we only converge the local potential and keep the non-local potential fixed to its last value. 

In the current implementation, you cannot specify when the inner loop is considered converged, we use a hard limit of 10'^-6^' with at most 50 iterations. The non-local potential is not mixed, we use only the new value. To switch to the nested self-consistency specify an \(imix \ge 10\). The last digit specifies the mechanism which is used for the convergence of the local potential.

Calculating with Hybrid Functionals
--------------
The input generator can be used to create a inp-file for hybrid functional calculations, too. For this you specify \(hybrid=T\) in the \(input\) namelist, which will generate a file inp_hyb that is used to run hybrid functional calculations. For further configuration options refer to the documentation of the input generator.

For the construction of the mixed basis, we need some starting wavefunctions. So you need to converge a calculation with any local functional before starting with hybrid functionals. Then you copy all necessary files (inp_hyb,enpara,sym.out,cdn1,kpts,eig,potcoul,pottot,potx,olap) to a new directory, where you rename the inp_hyb file for hybrid functional calculations to inp. Then start your calculation as usual.

When your calculation is finished, you can grep for the additional tag 'HF' in total energy and/or distance of charge to examine the convergence. Note that the total energy in the local iterations has no meaningful physical interpretation as the non-local part is missing. It is only printed to the output to examine the convergence of the local loops.


# InpHybFile

## Input file for hybrid functionals

    01|strho=F,film=F,dos=F,isec1=999,ndir= 0,secvar=F
     02|Silicon
     03|any any ,invs=T,zrfs=F,invs2=F,jspins=1,l_noco=F,l_J=F
     04|     0.00000000     5.13040000     5.13040000
     05|     5.13040000     0.00000000     5.13040000
     06|     5.13040000     5.13040000     0.00000000     0.00000000     1.00000000
     07|hse    non-relativi
     08|gcutm= 3.10000,mtol=0.00010000,lambda= 3,lexp=16,bands= 64
     09|
     10|**********************************
     11|Si  14    4    8  521  2.160000  0.022000
     12|
     13| 2,force =F,lcutm= 4,select= 4, 0; 4, 2,nlo= 0,llo=
     14|  1.000000  1.000000  1.000000  8.000000
     15| -1.000000 -1.000000 -1.000000  8.000000
     16|**********************************
     17| 11.000000  9.200000
     18|vchk=F,cdinf=F,pot8=T,gw=0,numbands= 70
     19|lpr=0,form66=F,l_f=F,eonly=F,eig66=F,soc66=T
     20|  6  8
     21|  1  0
     22|Window # 1
     23|  -2.80000  11.00000   8.00000
     24|   3.60000 =kmax
     25|gauss=F   0.00100tria=F
     26|  0.000000  0.000000,l_soc=F,spav=F,off=F
     27|frcor=T,slice=F,ctail=F,disp=F,kcrel=0,u2f=F,f2u=F,bmt=F
     28|itmax= 15,maxiter= 25,imix=17,alpha=  0.05,spinf=  2.00
     29|swsp=F  0.00
     30|lflip=F  1
     31|vacdos=F,layers= 0,integ=F,star=F,nstars= 0     0.00     0.00     0.00     0.00,nstm=0,tworkf=  0.000000
     32|iplot=F,score=F,plpot=F,band=F
     33|  0  0.000000  0.000000,nnne=  0,pallst=F
     34|xa=   2.00000,thetad= 330.00000,epsdisp=   0.00001,epsforce=   0.00001
     35|relax 111
     36|emin_dos=  -0.50000,emax_dos=   0.50000,sig_dos=   0.01500
     37|nkpt=   64,nx= 4,ny= 4,nz= 4,gamma=T
    

(7) Specifies the xc potential. You can choose between hse and pbe0. 

(8) Settings for the Hybrid functionals: 


|Key|Description
|--|--
| `gcutm`  | reciprocal cutoff for the interstitial plane waves of the mixed basis                                                                                                     |
| `mtol`   | tolerance value for the removal of linear dependencies from the MT part of the mixed basis                                         
| `lambda` | scaling parameter for the Ewald summation required in the calculation of the Coulomb matrix. A large (small) parameter lets the real (Fourier) space sum converge faster. |
| `lexp`   |  l cutoff for the expansion of plane waves inside MT spheres needed for the calculation of the Coulomb matrix                                                 

(13)   
`lcutm`: defines the angular momentum cutoff of the mixed basis for each atom type 

(18)   
`pot8=T` is required in a hybrid functional calculation to guarantee that the wavefunctions used in the first iteration are orthonormal. 

(20)   
For each atom type a l cutoff for the wave functions in the calculation of the wave function products is introduced. 

(27)   
Hybrid functionals for the core electrons are not implemented. If `fcore=T` the core states are not updated during the self-consistency. If `fcore=F` the core states are determined with respect to the PBE functional. 

(28)  
`imix`: For an favourable convergence of hybrid functional calculations, we introduce a *nested* convergence scheme, i.e. the Fock matrix is kept fixed until the density is converged with respect to the fixed Fock matrix. Then a new Fock matrix is setup etc. With imix=[10,13,15,17] you turn this nested convergence scheme on, where in the inner loops the density is mixed by straight/Broyden 1st/Broyden 2nd/ or Anderson mixing. 

(37)  
`gamma`: Specifies that the generated k-point set includes the \Gamma-point, which is required for hybrid functionals in order to treat the convergence of the Coulomb matrix at the \Gamma-point accurately. This option only works, if nx/ny/nz is set.

# Hands-on tutorial III: Approaches beyond LDA/GGA in the DFT self consistency

## Silicon

*approx. 60min* 

We already studied silicon intensively on the first day of our Hands-on seminar. Now, we focus on the changes which occur when using a hybrid functional for the calculation. In this example you will learn how the representation of the FLAPW basis functions via Kohn-Sham states converges with respect to the number of unoccupied bands used. If you have some spare time, you can investigate what convergence behavior the other new parameters show. 

**Goals:** 

*   generate input file for hybrid functional calculation for silicon 
*   converge a PBE calculation with a 4 \times 4 \times 4 k-point mesh 
*   examine the dependence of the bandgap on the number of unoccupied bands 

**Instructions** 

We want to examine a bulk sample of silicon with a rather small k-point mesh (to reduce the computation time). You can reuse the input file for the input generator from Monday. **Add the flag** `hybrid=t` **in the input namelist**. When you run the input generator you will notice that an additional file was generated the `inp_hyb` file. This can be used to calculate a hybrid functional calculation **starting from a converged calculation**. 

Let's have a closer look at the two input files. You will notice that the input file for PBE calculations differs slightly from the one you generated on Monday. The last line reads  
`nkpt=   64,nx= 4,ny= 4,nz= 4,gamma=T`  
now. This change is motivated by the fact that in hybrid functionals the \Gamma-point has to be included in the k-point mesh, because the Coulomb potential is divergent in its environment. The last switch ensures that this is always the case. **It is important that you calculate your PBE calculation with the same k-point mesh as your hybrid functional calculation**. 

The [`inp_hyb` file][2] shows some more changes. In the 7th line you observe that the name of the exchange-correlation functional is changed to `hse`. This starts calculations with the HSE06 functional. In this tutorial you will use the HSE functional or the PBE0 functional. In the latter case, you would replace the `hse` with `pbe0`. The next line gives you several parameters to converge the hybrid functional calculations. For this tutorial and also in general the number of `bands` is the most important one. It determines how many eigenstates are used to describe a FLAPW basis-function. There are some more changes in a few other lines, which we don't describe here, if you want to understand those changes please refer to [this page][2]. 

**Note:** There are two different parameters in the input file which determine the number of bands! The parameter `numbands` in line 18 determines how many eigenvalues and eigenvectors are written to the `eig` file. This parameter determines the maximum number of bands you can use for your calculation and it is present for PBE as well as hybrid functional calculations. The other parameter `bands` in line 8 determines the number of bands actually used for the transformation. It is only present for hybrid functionals and will determine whether your system is converged or not. This parameter has a significant impact on your computation time, so it is important to choose it appropriately. In the first part of this tutorial, we will study the effect of this parameter on the convergence of the bandgap. 

**Increase the number of bands (numbands, line 18) written to the output to 100**. Converge a calculation with PBE in the same way as in the first tutorial. If you compare your calculation to the one on Monday, you will notice that the `kpts` file is different, as the \Gamma-point is included as first k-point. There will be more eigenvalues written in the `out` file at every k-point. 

*We recommend to create a separate directory for every of the following calculations*. For a hybrid functional calculation you need all files generated by your PBE calculation. Additionally rename the `inp_hyb` file to `inp`. **Again, increase the number of bands (numbands, line 18) written to the output to 100**. For silicon 6 iterations should be sufficient. Choose one of the hybrid functionals (either PBE0 or HSE) and converge a silicon calculation using 20, 30, \ldots, 60 **bands** (line 8). Compare the direct bandgap at \Gamma as well as the indirect bandgaps \Gamma\rightarrow X and \Gamma \rightarrow L for the different numbers of unoccupied states. **Hint:** Silicon has four occupied bands and the internal coordinates of these three k-points are: 

| | | | |
|---|---|---|---|
| L      | 0.000 | 0.000 | 0.500 |
| \(\Gamma\) | 0.000 | 0.000 | 0.000 |
| X      | 0.000 | 0.500 | 0.500 |

To get the band transition, calculate the difference of the eigenvalue of the highest occupied state and the eigenvalue of the lowest unoccupied state at the respective k-points. 

Plot the different band transitions with a plotting program of your choice (gnuplot, xmgrace, ...). With how many bands are the band transitions converged? Compare the converged values with the result from PBE and experiment. 

**Advanced goals** (if you still have some time left) 

*   Is the result really converged, check using more bands. 
*   Calculate a convergence with respect to the parameter G_{cutm} - the planewave cutoff of the mixed basis. How much can you reduce the value, before the values of the band transitions change. 
*   Not all properties converge equally fast, consider for example the total energy, is it converged? Physically important are only total energy differences; try using a slightly different lattice constant and test the total energy difference between those two configurations is converged. 



* * *



## ScN

*approx. 60min* 

ScN crystallizes in rocksalt structure (fcc lattice with basis: one atom resides at the origin, where the other is shifted by half of the bodydiagonal) with a lattice constant of `8.50 a.u.` Experimentally, ScN is a semiconductor with an indirect band gap \( \left( \Gamma \rightarrow \mathrm{X} \right)\) of `1.3 eV`, and a direct gap \( \left( \Gamma \rightarrow \Gamma \right) \) of `3.8 eV`. 

*   Calculate the electronic ground-state with PBE, PBE0, and HSE! Compare the computed direct and indirect band gap for all three functionals! Choose a `4x4x4` **k**-point set and moderate parameters for the mixed product basis!  
    If the PBE calculation does not converge, remove the **Broyden** files! 
*   **Advanced:** Interpolate the PBE0 or HSE bandstructure of ScN along the **k**-point path  \mathrm{W}\rightarrow\mathrm{L}\rightarrow\Gamma \rightarrow \mathrm{X}\rightarrow \mathrm{W} \rightarrow \mathrm{K} employing the Wannier-interpolation technique!   
    In a first step try to interpolate the PBE bandstructure from the irreducible 4x4x4 **k**-point mesh. Check the quality of the interpolation by comparing with the PBE bandstructure that is obtained by exact diagonalization of the PBE Hamiltonian.  
    **Hint:** The **k**-points that define the path in the BZ are in internal coordinates given by: 

| | | | |
|---|---|---|---|
| W      | 0.250 | 0.500 | 0.750 |
| L      | 0.500 | 0.500 | 0.500 |
| \(\Gamma\) | 0.000 | 0.000 | 0.000 |
| X      | 0.000 | 0.500 | 0.500 |
| K      | 0.375 | 0.375 | 0.750 |  |

*[Step-by-Step-Guide][3] 

* * *



## Transition metal oxides

*approx. 120min* 

Localized states are very important in transition metal oxides because of the 3d-states. In the [LDA+U Method][4] one applies on-site Coulomb interaction to shift the states towards lower energies. However, typically the values for U are very different for the materials in the 3d series. In this part of the tutorial, you will investigate how the hybrid functionals with no tunable parameters perform for MnO. 

**Goals:** 

*   Calculate the direct band transition in MnO using the HSE functional (suggested k-point mesh 2\times2\times2). 
*   Determine the energy difference between the anti-ferromagnetic and the ferromagnetic state for an [111]- and an [001]-configuration. 
*   Calculate the J_{ij} values from these energy differences assuming a classical Heisenberg model. 
*   Compare the results with an LDA+U calculation. 

**Instructions** 

The ground state of MnO (a=4.445Å) is an anti-ferromagnetic configuration along the [111]-direction (ferromagnetic planes orthogonal to the space diagonal, adjacent planes couple anti-ferromagnetically) in the rocksalt structure. Setup the unit cell for the input generator and generate an input file for the PBE and HSE calculation. Remember to select an anti-ferromagnetic unit cell. [Image][6] [Hint][6] 

Remember to set the k-point mesh to \(2 \times 2 \times 2\) when you converge the PBE calculation and select then the HSE functional. Keep in mind that you should increase the number of bands written. Copy your converged PBE calculation in a new directory and converge the hybrid functional calculation (approx. 8-9 iterations). Check the band transition of your converged result. Do you find the insulating ground state? 

Now run a ferromagnetic calculation with the same unit cell for PBE and your selected hybrid functional. Check the energy difference between ferromagnetic and anti-ferromagnetic configuration. What is the ground state? 

Repeat analogously for an anti-ferromagnetism along \(\[001\]\) (in-plane ferromagnetic, out-of-plane anti-ferromagnetic) structure. Determine the energy difference in this configuration. Assume a classical Heisenberg model (i.e. |\mathbf{S}| = 1) 

\( H = -\frac 12\sum\_{i\neq j} J\_{ij} \mathbf{S}\_i \mathbf{S}\_j  \)

with only two different J's for nearest neighbors (J_1) and next nearest neighbors (J_2). Determine the coupling constants J_1 and J_2. [Hint][7] 

Repeat the anti-ferromagnetic calculation using a LDA+U calculation with the U and J values of [Ansimov *et al*.][9] (MnO U=6.9eV, J=0.86eV). Refer to the documentation of the LDA+U method for details. Compare the band transitions with the values of the hybrid functional. 

**Advanced goals** (if you still have some time left) 

*   Check the coupling constants with the LDA+U approach 
*   The hybrid functionals should be independent of the starting point. Confirm this by starting a hybrid functional calculation on top of your converged LDA+U calculation. Are the bandgaps and total energies the same? How many HF iterations does the calculation need to converge, how does this compare to a hybrid functional calculation starting from PBE? **Hint:** Remove the density matrix file before starting the hybrid functional calculation. 
*   Compare a PBE calculation of CoO with LDA+U (U=7.8eV J=0.92eV) and/or hybrid functionals. The t_{2g} states are only partially occupied and are degenerate in PBE. Is this degeneracy lifted in LDA+U and hybrid functionals? 



* * *



## InN

InN naturally crystallizes in the wurzite structure, but it can be synthesized in the zincblende structure, too. Here, we will focus on zincblende InN (fcc lattice with basis: one atom sits at the origin, the other one is shifted by a quarter of the bodydiagonal) with a lattice constant of `9.41 a.u.`.  
**Note:** For this example you have to compile the *Fleur* code without inversion symmetry and without SOC! 



*   Does PBE predict a semiconductor? 

*   Is InN a direct or indirect semiconductor in PBE0 and HSE? Choose a `4x4x4` **k**-point set!  
    **Hint:** If the hybrid functional calculation does not converge, describe the `4d` electrons of In with the LAPW basis and the `5d` ones with a local orbital. 

*   It is well known that LDA and PBE underbind `d` states with respect to experiment. The experimental energetic position of the `In` `4d`-states with respect to the Fermi level spreads from `-14.9 eV` to `-16.7 eV` in the literature. How does the energetic position of the occupied `d`-states of `In` change if going over from PBE to the hybrids? As an estimate evaluate the `d`-band position only at the \Gamma-point!

 []: http://www.flapw.de/pm/uploads/User-Documentation/afm1.pdf
 []: http://dx.doi.org/10.1103/PhysRevB.44.943
 
 
 # HybridGW

## Hybrid Functionals + GW

Algorithm to create the input for a GW calculation in a hybrid functional calculation 



1.  Converge a PBE (or any other local functional) calculation. The wave functions are necessary to create the mixed product basis as in normal hybrid functional calculations. 

2.  Run one iteration with this local functional and set `gw=1`. In this step, some lattice information is created, that does not depend on the functional. Furthermore the core states are stored, which will not be updated in hybrid functionals, because we employ the frozen core approximation for them. 

3.  Rename `inp_hyb` to `inp` and converge a hybrid functional calculation (`gw=0`). 

4.  Set `gw=2` and create the rest of the files needed for GW. **Please notice:** this step is not parallelized, yet, so you must run it on a single processor! 

5.  Proceed calculating the GW results with SPEX

 
 
