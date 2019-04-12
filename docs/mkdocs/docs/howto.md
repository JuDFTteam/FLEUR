# Relaxation

Here, we describe how to do a structural optimization. 

Suppose, you have converged a charge density, e.g. of a 3 layer Cu film as described in the example inp-file. Now, our interest is to optimize the position of the topmost layer (the position of the central layer is of course fixed by symmetry). What you have to do: 



## How to get forces 

*   change ` force =F ` to ` force =T ` in the inp-file, for those atoms, which you want to relax. 
*   change ` l_f=F ` to ` l_f=T ` in the inp-file, to allow the generation of Pulay-forces 
*   edit the line ` relax 000 001 ` at the end of the inp-file to allow relaxations in specific directions (here, ` 000` means no relaxation of the first atom, ` 001 ` means only in z-direction for the second). 
*   then run a few iterations, until the program stops with ` GEO: new inp created !`. This will happen, if the forces of two subsequent iterations do not differ more than 0.00001 htr/a.u. (This parameter should not be changed to ensure good convergence of the forces. It can, however, in cases of emergency be changed by creating an "eps_force"-file, which contains a different convergence parameter.) 

When this is finished, you will notice that you have two new files in your working directory: 

*   `inp_new` containing a new guess for the atomic positions and 
*   `forces.dat` which stores the positions and forces of this optimization step (and, if you keep it, also of the next steps). 
*   `out` contains the different Pulay force terms (cf. [![Symbol - externer Link][2]Rici Yu et al., PRB 43, 6411 (1991)][2] and the files `force*` of the Fleur source code). 



## How to start with a new geometry 

Now there are different ways how to proceed. In any case, you should not forget to remove a couple of files, that depend on the geometry: 

*   wkf2, the file containing the step-function (which cuts out the muffin-tins from the interstitial), and 
*   broyd and broyd.7, containing the (in/out)-density history of this self-consistency cycle. 

To be on the safe side, also the file "stars" can be removed, along with all "cdnXX", "pot[tot,coul]", "eig" and "tma[t,s]" files. Keep your directories clean. 

One possibility is to reuse just a minimum of files and create a new starting density and proceed as usual: 



*   ` mkdir ../GEO2 ` 
*   ` cp enpara kpts sym.out forces.dat ../GEO2 ` 
*   ` mv inp_new ../GEO2/inp ` 
*   ` cd ../GEO2 ` 
*   set ` strho=T ` 
*   run ` fleur.x ` etc. 

In other cases, it might be useful to reuse the old charge density with the new inp-file: 



*   ` mkdir ../GEO2 ` 
*   ` cp enpara kpts sym.out forces.dat cdn1 ../GEO2 ` 
*   ` mv inp_new ../GEO2/inp ` 
*   ` cd ../GEO2 ` 
*   ` echo " F" > qfix ` 
*   run ` fleur.x ` etc. 

In this case, you did not create a new starting density, but instead -- via the "F" in the file "qfix" indicate to the program that this is an old charge density, which is reused for shifted atomic positions (and therefore, the program renormalized the charge density a bit differently than normal). After the first iteration you will find that the "F" in this file changed to a "T". 



## Finishing the structural optimization 

This procedure can be repeated over and over again: 

1.  converge a calculation 
2.  obtain forces and an "inp_new" file 
3.  set up a new geometry and go to 1 

For once, you will notice that the "forces.dat" file gets longer and longer, containing the history of your structural optimization (like the broyd-files store your convergence of the electronic degrees of freedom). Again, a quasi-Newton scheme is used to minimize your forces. This has advantages, but can also lead to problems: 

*   if you are already in a parabolic region of your total-energy landscape, this scheme converges fast. Of course, this will depend on the number of degrees of freedom you have. 
*   if not, you might run into trouble. ` STOP 'bfgs: <p|y><0 ` is an indication of this case, but also very large relaxations (overlapping muffin-tins, atoms in the vacuum) can result. Then you should delete the "forces.dat" file. 

If you delete the "forces.dat" file after every step, you will effectively do a steepest-descent minimization. Since the relation between force and displacement is in general unknown, you can enter a Debye-temperature in the inp-file ("thetad= 330.000" near the end of the file). 

As in every good mixing method, you can also specify a mixing parameter for your quasi-Newton scheme (called "xa" in the inp-file, typically set to 2%). 

When you encounter a message "GEO: Des woas" (Austrian for: "that's it"), you either converged your system very well (either all forces are smaller than the parameter "epsforce" of the inp-file, or the displacements in the last step were smaller than "epsdisp), or all forces vanished for some other reason: Maybe, you set 

*   ` force =F ` for all atoms where forces can occur or 
*   ` relax 000 ` for directions where forces can be expected. 

In any case, you should check. If your forces are smaller than 0.001 htr/a.u., for most cases the structure can be considered relaxed.

 []: http://link.aps.org/doi/10.1103/PhysRevB.43.6411
 
 
