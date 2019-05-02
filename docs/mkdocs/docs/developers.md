
Contributors guide
======================================

Everyone is very welcome to contribute to the enhancement of FLEUR.
Please use the [gitlab service] (https://iffgit.fz-juelich.de/fleur/fleur) to obtain the
latest development version of FLEUR.


##Coding rules for FLEUR:
In no particular order we try to collect items to consider when writing code for FLEUR

- Instead of 'stop' use calls to judft_error, judft_warn, judft_end

- Do not read and write any files. Files are neither replacements for common-blocks nor a
storage for status variables.
 Only exceptions: 
+ you create nice IO subroutines in the subdirectory io
+ you write to the typical FLEUR output files

Useful info for developers
==============================================

## Using fleur with the HDF5 library and debugging it with valgrind

HDF5 has to be built with the same compiler that is also used to compile fleur. If adapted the following commands can be used to compile a HDF5 library for fleur:

+ 'FC=/usr/local/intel/impi/4.0.3.008/intel64/bin/mpiifort CC=/usr/local/intel/impi/4.0.3.008/intel64/bin/mpiicc CXX=/usr/local/intel/impi/4.0.3.008/intel64/bin/mpiicc ./configure --enable-fortran --enable-fortran2003 --enable-parallel --enable-using-memchecker --enable-clear-file-buffers'
+ 'make'
+ 'make install'
+ 'make check' (optional)

Note:

+ The paths have to be adjusted such that that compiler is used which is also used to compile fleur.
+ The parallel switch is not needed for every calculation: Only for parallel calculations in which HDF5 is also used for the eigenvector IO.
+ The last two command line switches in the configure command turn on initializations of irrelevant array parts in HDF5. If valgrind is not needed it is probably the better choice to leave them away. If left away valgrind will complain about missing initializations in the HDF5 library.
+ valgrind gives partially strange behavior if used together with the intel compiler. It would be better to use it together with gfortran.
+ At the moment HDF5 is needed in version 1.8.*. Usage of version 1.10.* yields some problems.

Furthermore to configure and start fleur with HDF5 the following has to be done:

+ In your .bashrc the HDF5 library has to be added to the LD_LIBRARY_PATH. This implies a line like 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/hdf5/current/hdf5/lib'
+ configure fleur with some line like 'CMAKE_Fortran_FLAGS="-I~/hdf5/current/hdf5/include" FLEUR_LIBRARIES="-L~/hdf5/current/hdf5/lib;-lhdf5_fortran;-lhdf5" ./fleur/configure.sh IFF'
