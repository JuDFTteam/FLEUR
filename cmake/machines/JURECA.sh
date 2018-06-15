if module list 2>&1 |grep -q -i intel
then
    echo "Intel toolchain used"
    FC=mpif90
    CC=mpicc
    FLEUR_LIBRARIES="-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64"
    
elif module list 2>&1 |grep -q PGI
then
    echo "PGI toolchain used"
    FC=mpif90
    CC=mpicc
    FLEUR_INCLUDEDIR="$XML2_ROOT/include"
    FLEUR_LIBRARIES="-lstdc++;-L$MKLROOT/lib/intel64;-lmkl_scalapack_lp64;-lmkl_intel_lp64;-lmkl_pgi_thread;-lmkl_core;-lmkl_blacs_intelmpi_lp64"
    
fi
