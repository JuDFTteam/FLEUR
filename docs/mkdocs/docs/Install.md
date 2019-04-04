Configuration and Installation of FLEUR
=========================================
We are aware of the fact that installing FLEUR can be a tricky task on many machines. While we tried to make the process
as userfriendly as possible, there are still a couple of challenges you might encounter. Please check with your system administrator to
see if all requirements for using FLEUR can be fulfilled on your system. For help register at the [MailingList](support.md) and post your questions there.
If you manage to compile on some system that can be of general interest, please consider to adjust the 'machines.md' file in the docs (Or report to fleur@fz-juelich.de if you do not know how to do that).

* [QuickInstall](#quick-guide)
* [Requirements](#requirements)
* [The configure.sh script & cmake](#configure)
* [How to adjust to your configuration](#how-to-adjust-the-configuration)
* [Running the automatic tests](#ci-tests)

#Quick guide 

If you are extremely lucky (and/or your system is directly supported by us) installation can be very simple:

* run the configuration script `'PATH_TO_SOURCE_CODE/configure.sh`. You can do that in any directory in which the 'build' directory should be created. The script accepts some useful arguments, you can run the script with `configure.sh -h`  to get a list of supported arguments.
* The script creates the build directory and runs cmake. If all goes well (look at the output) you can then change to the build directory and run `cd build; make`
* If make does not report any error you are done!

Please be aware that there are different executables that could be be build:

* `inpgen`: The input generator used to construct the full input file for FLEUR
* `fleur`: A serial version (i.e. no MPI distributed memory parallelism, multithreading might still be used)
* `fleur_MPI`: A parallel version of FLEUR able to run on multiple nodes using MPI.

Usually only the serial or the MPI version will be build. You can run the MPI-version in serial while it is of course not possible to use the non-MPI version with MPI.

You might want to [run the automatic tests](#ci-tests).

#Requirements 

There are a couple of external dependencies in the build process of FLEUR. 

**Required are:**

* *cmake*: The build process uses cmake to configure FLEUR. You should have at least version 3.0. Some options might require newer versions. Cmake is available for free at [www.cmake.org]([http://www.cmake.org).
* *Compilers*: You will need a Fortran compiler and a corresponding C-compiler (i.e. the two have to be able to work together via the iso-c bindings of Fortran). Please check our [list of compilers](#compilers) to see if your compiler should work.
* *BLAS/LAPACK*: These standard linear algebra libraries are required. You should try your best not to use a reference implementation from [Netlib](http://www.netlib.org) but look for an optimized version for your system. In general compiler and/or hardware vendors provide optimized libraries such as the MKL (Intel) or ESSL (IBM). If you do not have access to those, check [openBLAS]([http://www.openbas.net).
* *libxml2*: this is a standard XML-library that is available on most systems. If it is missing on your computer you should really complain with your admin. *Please note that you might need a development package of this library as well.* To compile this library yourself, see [xmlsoft.org](http://xmlsoft.org).

**Optional**:

FLEUR can benefit significantly if the following further software components are available. Please be aware that some of these can be difficult to use for FLEUR and see the [Instructions for adjusting your configuration](#configure) for details on how to provide input into the build process to use these libraries.

* *MPI*: Probably most important is the possibility to compile a version of FLEUR running on multiple nodes using MPI. If you have a proper MPI installation on your system this should be straightforward to use. 
* *HDF5*: FLEUR can use the HDF5 library for input/output. This is useful in two situations. On the one hand you might want to use HDF5 for reading/writing your charge density files to avoid having a machine-dependent format that can prevent portability. Also the HDF5 IO gives you some more features here. On the other hand you have to use parallel-HDF5 if you do IO of the eigenvectors in a MPI parallel calculation. This is needed if you can not store the data in memory or want to preprocess the eigenvectors. Please be aware that you need the Fortran-90 interface of the HDF5!
* *SCALAPACK/ELPA*: If you use MPI and want to solve a single eigenvalue problem with more than a single MPI-Task, you need a Library with a distributed memory eigensolver. Here you can use the SCALAPACK or [[http://elpa.mpcdf.mpg.de/|ELPA]] library. Please note that the ELPA library changed its API several times, hence you might see problems in compiling with it.
* *MAGMA*: FLEUR can also use the MAGMA library to run on GPUs. If you intend to use this feature, please get in contact with us.

You should also check the output of `configure.sh -h` for further dependencies and hints.

#Configure

The `configure.sh` script found in the main FLEUR source directory can (and should) be used to start the configuration of FLEUR. 
It is called as

 `configure.sh [-l LABEL ] [-d] [CONFIG]`

The most used options are:

* -l LABEL specifies a label for the build. This is used to custimize the build-directory to build.LABEL and can be used
to facilitate different builds from the same source.
* -d specifies a debugging build.
* CONFIG is a string to choose one of the preconfigured configurations. It can be useful if you find one which matches your setup.

More options are available. Please check the output of `configure.sh -h`

The `configure.sh` script performs the following steps:

1. It creates a subdirectory called 'build' or 'build.LABEL'. If this directory is already present, the old directory will be overwritten.
2. It copies the CONFIG dependent configuration into this directory (this is actually done in the script 'cmake/machines.sh'). The special choice of "AUTO" for CONFIG will not provide any further configuration but relies completely on cmake. You can specify a config.cmake file in the working directory (from which you call configure.sh) to modify this AUTO mode.
3 Finally cmake is run in the build directory to determine your configuration.


If you specify -d as argument of configure.sh, the string "debug" will be added to LABEL and a debugging version of FLEUR will be build, i.e. the corresponding compiler switches will be set.

You might want to check our page on 
[how to adjust to your configuration](InstallTroubles.md) if you run into trouble.

#How to adjust the Configuration
While `cmake` and the `configure.sh` script can determine the correct compilation switches automatically in some cases (mostly those known to us), in many other instances 
fine-tuning is needed. In particular you might have to:

* provide hints on which compiler to use
* provide hints on how to use libraries.

## Setting of the compiler to use
By defining the environment variables FC and CC to point to the Fortran and C compiler you can make sure that cmake uses the correct compilers. E.g. you might want to say

`export FC=mpif90`.

Please be aware that the use of CONFIG specific settings might overwrite the environment variable.

###Adding flags for the compiler
This should be done using the `-flag` option to `configure.sh`. So for example you might want to say `configure.sh -flag "-r8 -nowarn"`.

In general for a compiler [not known](#compilers) in cmake/compilerflags.cmake you need at least an option to specify the promotion of standard real variables to double precision (like the `-r8`). But additional switches can/should be used.

###Adding include directories
For libraries with a Fortran-90 interface, ELPA, HDF5, MAGMA, ... you probably will have to give an include path. This can
be achieved using the `-includedir` option. So you might want to say something like
`configure.sh -includedir SOMEPATH` 
###Adding linker flags
To add flags to the linker you can do

* add a directory in which the linker looks for libraries with `-libdir SOMEDIR`
* add the corresponding link option(s) with e.g. `-link "-lxml2;-llapack;-lblas"`. Please note that the specification is different from the compiler switches as different switches are separated by ';'.

### Further options:

There are more options you might find useful. Please
check `configure.sh -h` for a list.



##Compilers

FLEUR is known to work with the following compilers:

**INTEL**:

The Intel Fortran compilers (ifort) is able to compile FLEUR. Depending on the version you might experience the following problems:


1. Versions <13.0 will most probably not work correctly
2. Version 19.0 has issues with the debugging version of FLEUR.



**GFortran:**

GFortran is knwon to work with versions newer thant 6.3.


**PGI:**

The PGI compilers also can compile FLEUR. Here you need ad least version 18.4 but might still run into some problems.
  


#CI-Tests

After the build was finished you might want to run the automatic test. 

Just type `ctest` in the build directory for this purpose.

Please note:
* The tests run on the computer you just compiled on. Hence a cross-compiled executable will not run.
* You can use the environment variables `juDFT_MPI` to specify any additional command needed to start FLEUR_MPI. E.g. say `export juDFT_MPI="mpirun -n2 " `to run with
two MPI tasks.
* You can use the environment variable `juDFT` to give command line arguments to FLEUR. E.g. say `export juDFT='-mem'`.
* To run a specific test only (or a range of tests) use the `-I` option of ctest (check `ctest -h` for details)
* The tests are run in Testing/work. You can check this directory to see why a specific test fails.