# SsoHowto

## Spin-orbit coupling in spin spirals

The effect of spin-orbit coupling on spin spirals can be estimated by performing an ordinary spin-spiral calculation and afterwards applying spin-orbit coupling in a perturbative way, i.e. with Andersen's force theorem. 

In order to do this, follow the standard procedure described below or read [further options and explanations][1]. 



*   Note, that the combination of spin spirals and SOC does not work with several equivalent atoms per atom type, thus you might need to [break the symmetry][2] of the density first. 
*   Use a Fleur code compiled without `inversion`, with or without `soc`. 
*   In the '[nocoinp][3]' file set `l_ss=T`, `sso_opt=FFT`. 
*   In the '[inp][4]' file 
    *   Set `l_soc=T`. 
    *   In line 31, specify the angles theta and phi. These angles determine the rotation of the spin coordinate system (as used in the 'nocoinp' file) with respect to the real space coordinate system (as used to define the Bravais matrix). E.g. for theta=pi/2 and phi=0 you define a spin-rotation axis along x (in the 'inp' file the zenith angle theta is entered before the azimuth angle phi, whereas in the 'nocoinp' file the azimuth angle is entered before the zenith angle). 
    *   If needed, you can set `off=T` and explicitly specify which atoms are affected by SOC. 
*   Proceed as with a normal spin-spiral calculation (starting from a converged spiral density or by applying the force theorem; in the latter case you might like to set the magnetization in the interstitial and vacuum regions to zero by setting `<a class='wikilink' href='https://www.flapw.de/pm/index.php?n=User-Documentation.BmtHowto'>bmt</a>=T` in line 32 of the '[inp][4]' file). 
*   You can set `eig66=T` and use the eigenvector file(s) 'eig*' obtained by a previous spin-spiral calculation (and the density that was used to generate these eigenvectors). 
*   Using the 'eig*' files, Fleur calculates the spin-orbit expectation values of the occupied states and writes the corrected sum of eigenvalues in the file 'eigsso'. This file also provides the atom-resolved eigenvalue sums. 
    *   In this step, the program generates huge temporary files. Therefore, it might be necessary to enable the output buffer to save time (e.g. on the Jülich Jump system (IBM sp4) it can take hours(!) to write these files if buffering is disabled in $XLFRTEOPTS). 
    *   You can use k-point parallelization but cannot use eigenvalue parallelization. 

The method is explained in this [![Symbol - externer Link][6]paper][6].

 [1]: https://www.flapw.de/pm/index.php?n=User-Documentation.SsoDetails
 [2]: https://www.flapw.de/pm/index.php?n=User-Documentation.Symbreak
 [3]: https://www.flapw.de/pm/index.php?n=User-Documentation.Nocoinp
 [4]: https://www.flapw.de/pm/index.php?n=User-Documentation.Input
 []: http://dx.doi.org/10.1016/j.physb.2009.06.070
