
help=0
CLI_LIBDIR=""
CLI_INCLUDEDIR=""
while [ $# -gt 0 ]
do
  case "$1" in
    -h) help=1;;
	  -help|--help) help=1;; 
    -b) backup=1;;
    -backup) backup=1;;
    -g) gitupdate=1;;
    -c) shift;compiler=$1;;
    -dc) unset FC; unset CC; unset CXX ;;
  	-gitupdate) gitupdate=1;;
	  -l) shift;label=$1;;
	  -cmake) shift;cmake=$1;;
    -scalapack) shift;CLI_USE_SCALAPACK=$1;;
	  -hdf5) shift; CLI_USE_HDF5=$1;;
    -libxml2) shift; CLI_COMPILE_LIBXML=$1;;
	  -wannier) shift; CLI_USE_WANNIER=$1;;
    -edsolver) shift; CLI_USE_EDSOLVER=$1;;
	  -kplib)  CLI_USE_KPLIB=1;;
	  -mpi) shift; CLI_USE_MPI=$1;;
	  -magma) shift; CLI_USE_MAGMA=$1;;
	  -elsi) shift; CLI_USE_ELSI=$1;;
	  -gpu) shift; CLI_USE_GPU=$1;;
	  -chase) shift; CLI_USE_CHASE=$1;;
    -libxc) shift; CLI_USE_LIBXC=$1;;
	  -link) shift; CLI_LIBRARIES=$1;;
    -libdir) shift; CLI_LIBDIR="$CLI_LIBDIR $1";;
    -flags) shift; CLI_FLAGS=$1;;
    -includedir) shift; CLI_INCLUDEDIR="$CLI_INCLUDEDIR $1";;
    -elpa) shift;CLI_ELPA=$1;;
    -cmake_opts) shift;CMAKE_OPTIONS=$1;;
    -make) make_directly=1;;
    -ninja) use_ninja=1;;
    -warn_only) CLI_WARN_ONLY=1;;
    -d) debug=1;;
    -amd) CLI_PATCH_INTEL=1;;
    -spack) conf_spack=1;;
    -*) error="Unknown argument";;
    *)  break;;	# terminate while loop
  esac
  shift
done

if [ $# -gt 0 ]
then
    if [ "$label" = "" ]
    then
	label=$1;shift
    else
	error="You specified the -l switch and gave an additional LABEL argument"
    fi
fi
if [ $# -gt 0 ]
then
    error="Extra unknown arguments"
fi

#check if -h or  -help was given as argument
if [ $help -gt 0 ]
then
   echo -e "${RED}USAGE: configure.sh [options] [label]{$NC}"
   echo "
The following command line options are supported.
General options:
  -h            : print this help-page
  -c #          : ${GREEN}specify the compiler toolchain to use(details below)${NV}
  -dc           : delete FC,CC,CXX environment variables
  -l #          : label for the build. ${GREEN}It will be attached to the name of
                  the build directory.${NC}
                  This option can also be specified as last argument without -l
  -d            : build a debugging version of FLEUR (adds .debug to label)
  -g            : do a git pull first if this is a git version
  -t            : generate all tests including those that run longer
  -b            : backup an old build directory if present
  -make         : do not stop after configure script but run make directly
  -cmake #      : cmake executable to use
  -cmake_opts # : additional options for cmake can be specified here directly
  -amd          : apply some patches to the Intel MKL to run on AMD (very experiemental)
  -ninja        : use Ninja bild system instead of GNU make

Command line options to disable recommended libraries:
  -hdf5 false       : do not use HDF5. 
  -scalapack false  : do not use the SCALAPACK library

Command line options to switch on/off features. These options overwrite the results of
the test and might lead to the configuration to fail.
  -elsi     [TRUE|FALSE} : use the ELSI library
  -wannier  [TRUE|FALSE] : use Wannier90 library
  -mpi      [TRUE|FALSE] : compile the MPI parallel version
  -libxc    [TRUE|FALSE] : use libxc library
  -edsolver [TRUE|FALSE] : use the Exact Diagonalization library by Jindrich Kolorenc
  -libxml2   true        : try to download libxml2 and compile it (experimental)
  -magma     true        : use the Magma library (no test,experimental)


Command line option to compile external libraries:
  -external # : download and compile external libraries before building FLEUR
                currently 'xml2', 'elpa' and 'chase' are possible options. The switch
                can be specified multiple times

Options to specify Fortran/Linker flags. Usually it is better to use the enviroment variables as
given below.
  -link #       : String to use for linking (options separated by ;, e.g. '-lxml2;-lhdf5')
  -libdir #     : Directory to find libraries in (can be specified multiple times)
  -flags #      : String to add while compiling (e.g. '-g')
  -includedir # : Directory to find include files (can be specified multiple times)

Special options:
  -gpu # : Compile for GPU. Currently you should specify something like acc:cc80 to use OpenACC
           and NVIDIA compute capability 80. Currently this works only using the NVHPC compilers.


To help the script finding a proper configuration you should provide some information on your compiler toolchain.
In detail, the script will determine the Fortran, C and C++ compilers to use. Best practise is to set the following
environment variables:
  FC         -- name of Fortran compiler
  CC         -- name of C compiler
  CXX        -- name of C++ compiler

Alternatively, you can try to specify the compiler toolchain using the -c option. The following options are available:
  mpi        -- use standard mpi-wrapper scripts (mpifc,mpicc,mpicxx). 
                This is a good option if you have a correct MPI installation available.
  intel      -- use the Intel oneAPI compilers (ifx,icx,icpx).
  intel-old  -- use the classic intel compilers (ifort, icc, icc)
  nvidia     -- use the NVHPC compilers
  gfortran   -- use the compilers from the GNU collection

  auto       -- let cmake determine the compilers (will usually use the 'first' compiler found if multiple are available)
  interactive-- search for compilers on you machine and offer a choice. This will also be offered if you specify no
                compiler toolchain and you have not set all the environment variables (FC,CC,CXX)
               
"
  exit 1
fi

if [ "$error" != "" ]
then
    echo $error
    echo "${RED}ERROR in calling configure. Please use -h for help on usage${NC}"
    exit 1
fi
