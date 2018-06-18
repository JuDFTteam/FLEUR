if module list 2>&1 |grep -q -i CMake
then
if module list 2>&1 |grep -q -i intel
then
    echo "Intel toolchain used"
    export FC=${FC:=mpif90}
    export CC=${CC:=mpicc}
    #determine XML2 module
    xml2=`module --show_hidden spider libxml2 2>&1 |grep libxml2/|grep -v module |tail -1`
    module load $xml2
    
    #determine ELPA module
    elpa=`module spider ELPA 2>&1 |grep hybrid`
    module load $elpa
    CLI_ELPA_OPENMP=1
    FLEUR_LIBDIR="$FLEUR_LIB $ELPA_LIB"
    FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $ELPA_MODULES_OPENMP"
    #load hdf5 module
    module load HDF5
    FLEUR_LIBDIR="$FLEUR_LIBDIR $HDF5_DIR/lib"
    FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $HDF5_DIR/include"
 
 
    FLEUR_LIBRARIES="-lelpa_openmp;-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64"
    ml 
elif module list 2>&1 |grep -q PGI
then
    echo "PGI toolchain used"
    FC=mpif90
    CC=mpicc
    FLEUR_INCLUDEDIR="$XML2_ROOT/include"
    FLEUR_LIBRARIES="-lstdc++;-L$MKLROOT/lib/intel64;-lmkl_scalapack_lp64;-lmkl_intel_lp64;-lmkl_pgi_thread;-lmkl_core;-lmkl_blacs_intelmpi_lp64"
else
    echo "You need to load the modules for the compiler"
    echo "e.g. module load intel-para"
    exit
fi
else
    echo "You should load the CMake module and the modules for the compiler"
    echo "e.g. do a 'ml intel-para CMake'"
    exit
fi
