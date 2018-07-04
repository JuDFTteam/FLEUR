#Configuration for JURECA@FZJ; use PGI toolchain
#determine XML2 module
xml2=`module --show_hidden spider libxml2 2>&1 |grep libxml2/|grep -v module |tail -1`
#determine ELPA module
elpa=`module spider ELPA 2>&1 |grep hybrid`

if module list 2>&1 |grep -q PGI
then
    echo "PGI toolchain used"
    FC=mpif90
    CC=mpicc
    FLEUR_INCLUDEDIR="$XML2_ROOT/include"
    FLEUR_LIBRARIES="-lstdc++;-L$MKLROOT/lib/intel64;-lmkl_scalapack_lp64;-lmkl_intel_lp64;-lmkl_pgi_thread;-lmkl_core;-lmkl_blacs_intelmpi_lp64"
else
    echo "You need to load the modules for the compiler"
    echo "e.g. module load intel-para CMake HDF5 $xml2 $elpa "
    exit
fi
