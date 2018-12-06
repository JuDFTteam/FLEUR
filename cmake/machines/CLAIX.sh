#Configuration for CLAIX@RWTH; do a 'module load intelmpi' before

  if  module list 2>&1 | grep -q intel
  then
    if  module list 2>&1| grep -q openmpi
    then
      echo "Please use intelmpi, e.g. do a module switch openmpi intelmpi"
      exit	
    fi
    if module list 2>&1| grep -q intelmpi
    then
       export FC=$MPIFC
       export CC=$MPICC
    fi
  elif  module list 2>&1 | grep -q pgi
  then
    echo "found PGI compiler"
  fi

  #ELPA
  if [ $ELPA_MODULES ]
  then
    CLI_ELPA_OPENMP=1
    FLEUR_LIBDIR="$FLEUR_LIBDIR $ELPA_LIB"
    FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $ELPA_MODULES"
  fi

 
