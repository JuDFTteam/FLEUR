Guide to Local Orbitals (v26)
==========================

## When do I need LO's? 

When the accurate description of an atom requires two 'valence' bands with equal symmetry (e.g. 3p and 4p in Ti compounds) then the use of a 'second window' (energy panel) or of LO's is recommended. There are three categories of problems that might require additional local orbitals: 


### (a) Atoms that tend to give raise to 'ghost-band' problems 

Some materials with high-lying core states (say, 5p in Cs) tend to produce ghost-band, i.e. the unwanted appearance of these high-lying states in the valence window. Also in more robust atoms, the use of small muffin-tin radii (that makes the basis more flexible) can give such problems. A typical precursor of a ghost-band is the output 'n eigenvalue(s) below the energy e\_low', with n >= 1 and a reasonably chosen lower energy bound e\_low. When the ghost-band appears in the valence window, i.e. its energy is higher than e\_low, the charge-density distance 'jumps' to a large value and one of the energy parameters drops to a value near e\_low. The band associated with this energy parameters should then be described with LO's. 



### (b) Low-lying bonding orbitals

Sometimes an atom binds with low- and high lying states of equal symmetry (e.g. in Ti oxides, the 3p binds to the oxygen while the 4p is involved in metallic states). Then one of these states should be described with LO's. 


### (c) Extreme spin-splitting of f-bands 

An extremely spin-split valence band, e.g. the 4f of Gd, requires two energy parameters to be described correctly. Here, LO's can be used for one spin channel.

 
## How do I calculate with LO's? 

Having identified the atom and state that should be calculated with a LO, two steps are necessary to switch from a normal to a LO calculation: 


### (a) Specify the number and l-character of the LO's (nlo and llo) and subtract them from the core-levels and adjust the number of electrons in the valence band. 

First, Identify the atom in the input file. Say, an iron has a very small muffin-tin radius and gives 3p ghost-bands: 


```
    | *************************                  |
    | Fe  26    7    8  461  1.700000  0.023000  |
    |                                            |
    |  1,force =F,nlo= 0,llo=                    |
```    


*   Now, we want a single LO: set nlo= 1 
*   it should be a p-orbital: set llo= 1 (use the l-quantum number) 
*   subtract two core-levels: 7 -> 5 (the 3p1/2 and 3p3/2) 

This gives us: 


```
    | *************************                  |
    | Fe  26    5    8  461  1.700000  0.023000  |
    |                                            |
    |  1,force =F,nlo= 1,llo= 1                  |
```    

Finally, add the electrons to the valence band, in our example change from: 


```
    | Window # 1                     |
    |   -0.50000   0.80000   8.00000
```    

to: 


```
    | Window # 1                     |
    |   -2.00000   0.80000  14.00000 |
```    
(the 3p contains 6 electrons, therefore 8+6=14. Here, we also had to
change e_low from -0.5 to -2.0, so that the 3p fits in the valence
window.)
    

### (b) Choose the energy parameter and fix these values or skip the states in the calculation of the energy parameters. 

Now set the energy parameters in enpara. (If you start a calculation, this should be done automatically, but check, whether it did, what you intended!) Taking the Fe-example from above, your enpara file probably contained energy parameters like: 


```
    |      energy parameters for window 1 spin 1                            |
    |      atom     s        p        d        f                            |
    |  --> 1    0.13301  0.21323  0.24802  0.24094 change: TTTT skiplo:   0 |
```    

and you expect the 3p somewhere near -1.7 htr., so you add: 


```
    |  --> 1    0.13301  0.21323  0.24802  0.24094 change: TTTT skiplo:   3 |
    |  --> lo  -1.70000                                                     |
    |  --> change   T                                                       |
```    

Here, in the last line we specified 'change T ', i.e this LO energy parameter should be adjusted to the center of gravity of the p-band. But, since we have two p-bands, we also want two energy parameters ! Therefore we specified in the first line 'skiplo: 3', that means: skip the first three bands for the evaluation of the 4p energy parameter. Alternatively, we could have set 'change: TFTT skiplo: 0', then the first p-energy parameter (0.21323) would have kept fixed at this value. In general, you have to avoid two equal energy parameters. Otherwise, the basis gets overcomplete and the program crashes. 



### Summary: 

*   set number and type of LO's 
*   reduce the number of core levels 
*   add valence electrons 
*   adjust lower energy bound of window 
*   specify energy parameter 
*   choose, which energy parameters to change 
*   set skiplo

# Restrictions and problems 

Note, that using the LO's is incompatible with the following features of the program: 



*   second-variation and Wu-diagonalization 
*   constrained magnetic moments calculations 

Sometimes, we can solve our problems by making a 2-window calculation rather than a local-orbital one. We are presently working on a version that includes (a) and (b) with LO's. We should also note that no explicit use of z-reflection symmetry in the diagonalization is made when LO are specified. 



## Common problems: 

*   After a calculation without LO's, in the first iteration with LO's, the charge density distance 'jumps' to a very high value. 
    
    *   Either the number of core levels was not changed correctly or the energy window was too small to include the LO's. 
    *   In the latter case "n eigenvalue(s) below the energy e\_low" can be found in the output-file: change e\_low. If that is not the case, check the number of core-levels and valence electrons. 
    

*   Program stops with 'rdc_22-notpos' or 'franza1' 
    
    *   Probably two energy parameters of equal 'l' converged to the same value (e.g. 3p and 4p). -> Fix one of them or set skiplo to correct value. 
    

*   The calculation does not converge. 
    
    *   Sometimes, the "normal" basisfunctions and the LO's cannot decide, who describes upper and lower band. Then the energy parameters of the normal basis and LO's alternate. 
    *   Fix both energy parameters to reasonable values. Finally, you can try to set 'change' to 'T' again.

# APW+LO
The APW+LO method was introduced by Sjöstedt et al. [Solid State Comm. 114, 15 (2000)](http://dx.doi.org/10.1016/S0038-10989900577-3). Formally, this new basis set can be constructed by dropping the "B" coefficients of the LAPW inside the sphere (thus creating a 'kink' of the wavefunction at the muffin-tin boundary) and linearising this APW-like basis set by a local orbital constructed from the energy-derivative of the radial basisfunction (the u-dot that appears in the LAPW basis set associated with the "B" coefficient). Thus, three things have to be done: 


*   (a) replace the kinetic energy operator (as in APW) 
*   (b) remove the "B" coefficients 
*   (c) add (2l+1) local orbital for relevant l's 

(a) is done by compiling the code with CPP_APW. The use of the APW-like kinetic energy operator is completely general and can also be used for LAPW calculations. 

(b) is done automatically for a atom where a 'u-dot'-LO is introduced. In the following we call these local orbitals "lo's", to distinguish from the "normal" LO's. 

(c) is done by specifying "negative" l's for the lo's, e.g. -1 for s, -2 for p etc.; the energy-parameter has to be given in the "enpara" file but is set to the same value as otherwise they are treated as LO's. 

Example 1: To calculate Ti with the APW+LO method, we specify 3 lo's (s,p,d) in the inp-file: 


```
    | *************************                  |
    | Ti  22    4    8  361  2.500      .03      |
    | ...                                        |
    |  1,force =F,nlo= 3,llo=-1-2-3              |
```    

and the enpara file is accordingly: 


```
    energy parameters for window 1 spin 1
         atom     s        p        d        f
     --> 1   -0.27040 -0.22782 -0.20974 -0.20411 change: TTTT skiplo:   0
     --> lo  -0.27040 -0.22782 -0.20974
     --> change   T        T        T    
```    

Now, this input can be used with a programm compiles with -DCPP_APW. 

Example 2: If you want to use additional LO's (in the above example for the 3s and 3p states of Ti), the input reads: 


```
    |  1,force =F,nlo= 5,ll0=-1 0-2 1-3        |
```    

and (supposing the 3s is at -2.3 htr and the 3p at -1.4): 


```
    --> 1   -0.27040 -0.22782 -0.20974 -0.20411 change: TTTT skiplo:   4
     --> lo  -0.27040 -2.30000 -0.22782 -1.40000 -0.20974
     --> change   T        T        T        T        T
```    

Example 3: For the accurate calculation of the forces an additional "polarization orbital", in the case of Ti this would be a f-orbital, is necessary. Instead of adding another lo, you can use a f-LAPW by changing the line containing the l-cutoffs for the Hamiltonian, e.g. from 


```
    |    8                                     |
```    

to 


```
    |   38                                     |
```    

where the additional '3' indicates that every orbital up to l=3 that is not treated as lo should be an LAPW. In our case this yields: 


```
    l      type
    
     0      APW+lo
     1      APW+lo
     2      APW+lo
     3      LAPW
     4      APW
     ...
     8      APW
```    

This implementation was done in collaboration with G. K. H. Madsen. More details and applications can be found in his paper [Phys. Rev. B **64**, 195134 (2001)](http://10.0.4.79/PhysRevB.64.195134).

 
