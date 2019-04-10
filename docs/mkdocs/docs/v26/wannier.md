# Installation

<< [Wannier Functions][1] | > [Step-by-step guide][2] 



## Installation

*   See the [download][3] page for a version of Fleur including Wannier functions. 

*   The Library of the Wannier90 Program `libwannier.a` has to be linked to the Fleur program during compilation. Download the source of Wannier90 (version v1.1) at [![Symbol - externer Link][5]http://www.wannier.org/][5] and compile both the Wannier90 program and the library `libwannier.a`. As described in the wannier90 documentation, you need a file `make.sys` corresponding to the compiler and libraries specific to your computer. However, if you use the Juelich computer cluster or Juelich super computers, you may also find appropriate files `make.sys` [here][5]. 

*   Add the line "#define wann" into the file "Imakefile" and run `imake`. Then run `make` to compile `FLEUR`. 

*   For the calculation of Wannier functions, the full Brillouin zone is needed. To generate the corresponding k-point files, it is sometimes convenient to use an external k-point generator instead of the one included in the Fleur program. You may generate a k-point set using the following Fortran programs: [kpointgen][6] generates the kpts file and [w90kpointgen][7] generates the [w90kpts][8] file.

 [1]: https://www.flapw.de/pm/index.php?n=User-Documentation.WannierFunctions
 [2]: https://www.flapw.de/pm/index.php?n=User-Documentation.Step-by-stepGuide
 [3]: https://www.flapw.de/pm/index.php?n=FLEUR.Downloads
 []: http://www.wannier.org/
 [5]: https://www.flapw.de/pm/index.php?n=User-Documentation.Wannier90-installation-make-sys
 [6]: https://www.flapw.de/pm/index.php?n=User-Documentation.Wannier-kpointgen
 [7]: https://www.flapw.de/pm/index.php?n=User-Documentation.Wannier-w90kpointgen
 [8]: https://www.flapw.de/pm/index.php?n=User-Documentation.Wannier-w90kpts
 
 
