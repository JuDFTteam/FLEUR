# LDAu

To introduce a Hubbard U for the description of strongly correlated electrons (e.g. in transition-metal oxides or f-electron systems) it is possible to use the LDA+U method (see e.g. [Anisimov et al. J.Phys.: Condens. Matter 9 (1997) p. 767-808][2]). It was implemented similar to [Shick et al., Phys. Rev. B 60, 10763 (1999)][3] but without the pseudo-perturbative treatment described there. For the input of the U and J parameters for an atom the second line (previously containing the local symmetry) is used: 



    ***************************
     Gd  64 ...
     &ldaU l=3,u=6.7,j=0.7 /
    

This is a namelist input, the format is free. The 3 parameters, l, u and j specify the orbital (l=0 for s, 1 for p, 2 for d and 3 for f), the Hubbard U and exchange interaction parameter J (given in eV). The latter two parameters can be taken from literature or be calculated [Solovyev and Dederichs, Phys. Rev. B 49, 6736 (1994)][4]. The density matrix is stored on a file [n\_mmp\_mat] and should be kept with the charge density.  
For vi-users, a convenient way to add LDA+U to all atoms of a specific type is (in the above example) 

    :g/Gd  64/+1 s/ /\&ldaU l=3,u=6.7,j=0.7 \//
    

Since the double-counting correction of the LDA+U method can be formulated in different ways, the atomic and the "around mean-field" limes, a switch has been introduced in the line 



    &ldaU l=2,u=4.0,j=2.0,l_amf=T/
    

If l_amf=T is set, the around mean-field limes is used, otherwise the atomic limes is taken. For practical purposes the difference is, that in the atomic limes (default in the version 0.22) states that are more than half occupied are lowered in energy, while in the around mean-field limes states that are more occupied than the average are lowered in energy. For a comparson of the two methods and another variant (not included) see e.g. [Petukhov et al., PRB 67, 153106 (2003)][6]. 

**Using LDA+U:** 

(1) *Setting up the density matrix:* When you have converged a charge density with LDA or GGA, run a single iteration with `&ldaU ...` in the inp-file. A line in the output will inform you that LDA+U has been skipped, since no density-matrix was found, but together with the charge-density a `n_mmp_mat` file will be created, that contains the density-matrix to start with. Then delete the files `broyd` and `broyd.*`. 

(2) *Choosing the right executable:* Note, that if the density matrix contains off-diagonal elements you will get a complex-hermitian Hamiltonian, e.g. if a orbital moment arises. Then, the program should be compiled without inversion, even if inversion-symmetry is present in the lattice. 

(3) *Converging the calculation:* If you now run the executable, the charge-density and density-matrix have to be converged again. By default, the density-matrix is converged separately using straight-mixing with the mixing parameter given in the last line of the `n_mmp_mat`-file. If you delete this line, the density-matrix is mixed together with the charge-density in the `broyd`-files. If you switch from straight- to broyden-mixing, delete these files. 

(4) *Finishing:* Monitor the convergence of the charge-density and density-matrix by e.g. `grep dist out`. Note that both densities have to be converged, which is sometimes not trivial. Also note, that the final result may depend on the starting point, e.g. different orbital momenta can be stabilized if you calculate an isolated atom starting from different density-matrices. 

**Limitations and Problems** 



*   Only one set of orbitals (one `l`-value) can be specified per atom (anyway, you should use LDA+U with care and only for those cases, where it is really needed). 
*   If local orbitals are present with the same `l`-value as the one where the U acts on, the method will not discriminate between the local orbital and the normal LAPW's. 
*   Forces calculated for an atom, where the LDA+U method is applied, are currently inaccurate (of course this depends on the variation of the occupation matrix with the atomic displacement). See the discussion in [Tran et al., Comp. Phys. Comm. **179**, 784 (2008)][7]. 
*   In principle, there exists also a spin-offdiagonal part of the density-matrix, which is not implemented so far (P. Novak, A. Shick et al.)

 [2]: http://dx.doi.org/10.1088/0953-8984/9/4/002
 [3]: http://link.aps.org/abstract/PRB/v60/p10763
 [4]: http://link.aps.org/abstract/PRB/v49/p6736
 [6]: http://link.aps.org/abstract/PRB/v67/e153106
 [7]: http://dx.doi.org/10.1016/j.cpc.2008.06.015
