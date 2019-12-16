
  if module list 2>&1| grep -q mpi.intel
  then
    export FC=mpif90
    export CC=mpicc
  fi

  if  module list 2>&1 | grep -q elpa
  then
    #ELPA
    CLI_ELPA_OPENMP=1
    FLEUR_LIBDIR="$FLEUR_LIBDIR $ELPA_LIB"
    FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $ELPA_MODULES_OPENMP"
  fi

