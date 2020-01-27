
  if  module list 2>&1 | grep -q intel
  then
    if  module list 2>&1| grep -q mpi.ibm
    then
      echo "Please use intelmpi, e.g. do a module switch mpi.ibm mpi.intel"
      exit	
    fi
    if module list 2>&1| grep -q mpi.intel/2017
    then
      echo "There are problems with mpi.intel/2017, try another version"
      exit
    fi
    if module list 2>&1| grep -q mpi.intel
    then
      export FC=mpif90
      export CC=mpicc
    fi
  fi

  if  module list 2>&1 | grep -q elpa
  then
    #ELPA
    CLI_ELPA_OPENMP=1
    FLEUR_LIBDIR="$FLEUR_LIBDIR $ELPA_LIB"
    FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $ELPA_MODULES_OPENMP"
  fi

