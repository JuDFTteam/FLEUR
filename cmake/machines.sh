#This file contains the defaults for compiling FLEUR on known machines

#please add further definitions here and also add code below!
read -r -d '' known_machines << EOM
   IFF         -- cluster @ PGI of FZJ
   JURECA      -- @JSC
   JUQUEEN     -- @JSC
   CLAIX       -- @RWTH
   MARCONI     -- @CINECA 
EOM


function configure_machine(){
    
    if [ "$machine" = "JURECA" ] 
    then
	echo "JURECA configuration used"
	if module list 2>&1 |grep -q intel-para
	then
	    echo "Intel toolchain used"
	    if module list 2>&1| grep -q Python/2.7.12
	    then
		echo "Python module loaded for XML (OK)"
	    else
		echo "You have to load the Python module"
		echo "module load Python/2.7.12"
		exit
	    fi
	    module load ELPA/2016.05.003-hybrid
	    module load HDF5
	    export CC=mpicc
	    export FC=mpif90
	    export CMAKE_Fortran_FLAGS="$CMAKE_Fortran_FLAGS -I$EBROOTELPA/include/elpa_openmp-2015.11.001/modules -I$EBROOTHDF5/include -mkl"
	    export FLEUR_LIBRARIES="$FLEUR_LIBRARIES;-L$EBROOTELPA/lib;-lelpa_openmp;-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64;-L$EBROOTHDF5/lib;-lhdf5;-lhdf5_fortran"
	elif module list 2>&1 |grep -q PGI
	then
	    echo "PGI toolchain used"
	    echo "Needs to be configured"
	    exit
	else
	    echo "You have to load the correct modules for compiling"
	    echo " a) intel-para, python/1.7.12"
	    echo " or"
	    echo " b) PGI"
	    exit
	fi
    # JUQUEEN
    elif [ "$machine" = "JUQUEEN" ]
    then 
   	echo "JUQUEEN configuration used"
	module load  hdf5/1.8.15_BGQ
	module load  scalapack/2.0.2_elpa_simd
	export CC=mpixlc
	export FC=mpixlf2008_r
	export CMAKE_Fortran_FLAGS="$CMAKE_Fortran_FLAGS -I${HDF5_DIR}/include -I${ELPA_INCLUDE}"
	export FLEUR_LIBRARIES="$FLEUR_LIBRARIES;-L$SCALAPACK_ROOT/lib;-lelpa;-lscalapack;-L/bgsys/local/lapack/3.3.0_g/lib;-llapack;-L/bgsys/local/lib;-qessl;-lesslsmpbg;-L$XML2LIB;-lxml2;-L${HDF5_DIR}/lib;-lhdf5_fortran;-lhdf5;-L/bgsys/local/zlib/lib/;-lz;-L/bgsys/local/szip/lib/;-lsz"
	
    #IFF linux cluster
    elif [ "$machine" = "IFF" ]
    then
	echo "IFF cluster configuration used"
	export FC=mpiifort
	export FLEUR_LIBRARIES="$FLEUR_LIBRARIES;-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64"
    #RWTH cluster
    elif [ "$machine" = "CLAIX" ]
    then
	echo "CLAIX@RWTH configuration used"
	if ! module list 2>&1| grep -q intelmpi
	then
	    echo "Please use intelmpi, e.g. do a module switch openmpi intelmpi"
	    exit
	fi
	module load LIBRARIES
	module load hdf5
	export FC=mpiifort
	export FLEUR_LIBRARIES="$FLEUR_LIBRARIES;-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64"
    elif [ "$machine" = "MARCONI" ]
    then
	if ! module list 2>&1| grep -q " intel\/" || ! module list 2>&1| grep -q " intelmpi" ||! module list 2>&1| grep -q cmake
	then
	    echo "Load the modules needed to compile: intel,intelmpi,cmake"
	    exit
	fi
	export FC=mpif90
	export FLEUR_LIBRARIES="$FLEUR_LIBRARIES;-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64"
    elif [ "$machine" = "AUTO" ] 
    then
	echo "No machine specific settings used"
	echo "GOOD LUCK!"
    else
	echo "No valid machine configuration specified"
	exit
    fi
}

