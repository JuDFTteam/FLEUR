# Magnetism

## starting a magnetic calculation 

When you use the input-generator, this program tries to set up a magnetic calculation, whenever atoms like Fe, Co, or Ni appear in the input. Then ` jspins=2 ` is set in the inp-file and at the bottom a line like 

    swsp=F  2.20  1.60  0.90
    

appears, indicating that the starting density should have magnetic moments of 2.2, 1.6 and 0.9 Bohr for the atom types 1, 2 and 3. If you want to start with different moments, e.g. an antiferromagnetic calculation, change this line like 

    swsp=F  2.20 -2.20  0.00
    

then atoms 1 and 2 have antiparallel spins. Of course, the input generator does not know, that NiAl is nonmagnetic, Gd2Fe is antiferromagnetic, or a thin film of Rh is magnetic. This you have to specify yourself. 



## spin-polarizing a calculation 

If you started from a non-magnetic calculation and want to spin-polarize your calculation, 



*   make sure that you have a file "cdnc" in your working directory and set ` jspd=2 ` in the fl7para file (or simply remove it), 

*   set `jspins=2` and ` swsp=T ` in the inp-file, followed by the magnetic moments you want to induce in the atoms (see above). 

*   running fleur.x results in STOP 'spin polarized density generated' 

*   Now set again ` swsp=F ` in the inp-file and converge 

If you want to flip the magnetization of some atom n in a spin-polarized calculation, set l\_flip=T and set the n'th digit in the same line to -1. Run fleur.x and after the STOP set l\_flip=F again.

# Non-CollinearMagnetism

A short primer in how to do calculations using Non-Collinear Magnetism. 



1.  First of all, converge a collinear ferromagnetic configuration. 
2.  Then you should copy cdn1 to rhomat_inp. 
3.  Next create a nocoinp-file (make sure, that it is in the right format!). If you want to use the constraint it is probably best to switch it on right from the start (l\_mperp=T,l\_constr=T). A good choice for mix_b is 0.5. 
4.  You should now delete fl7para and the broyd-files. 
5.  Then set l_noco=T in inp. 
6.  If ctail=T in inp, set ctail=F. 
7.  Then converge as usual noticing that rhomat_inp is now used instead of cdn1. 

** Be careful** 

when you set `l_ss=T` ! 

1.  This cannot be combined with constraints (`l_constr=T`) because it's not implemeted yet. 
2.  When setting `l_soc=T` a perturbative scheme is used as the generalized Bloch theorem does not hold in the presence of spin-orbit coupling. 

when you set `l_mperp=T` ! Then two things do not work correctly: 

1.  the orbital moments in the out- and inf-file get wrong 
2.  the density matrix (in file `n_mmp_mat`) is wrong 

these things still have to be implemented to work together. 

For an introduction into the technique, see the Psi-k Highlight [FLAPW goes Non-collinear][4]

 [4]: http://psi-k.dl.ac.uk/newsletters/News_38/Highlight_38.pdf
