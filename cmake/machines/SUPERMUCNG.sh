
  if module list 2>&1| grep -q mpi.intel
  then
    export FC=mpif90
    export CC=mpicc
    export CXX=mpicxx
  fi

  if  module list 2>&1 | grep -q elpa
  then
    #ELPA
    CLI_ELPA_OPENMP=1
    FLEUR_LIBDIR="$FLEUR_LIBDIR $ELPA_LIB"
    FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $ELPA_MODULES_OPENMP"
  fi

  if  module list 2>&1 | grep -q hdf
  then
    CLI_USE_HDF5=$1
    FLEUR_LIBDIR="$FLEUR_LIBDIR $HDF5_BASE/lib"
    FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $HDF5_BASE/include"
  fi

