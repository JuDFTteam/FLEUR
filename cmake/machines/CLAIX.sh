    if ! module list 2>&1| grep -q intelmpi
    then
      echo "Please use intelmpi, e.g. do a module switch openmpi intelmpi"
      exit	
    fi


    export FC=$MPIFC
    export CC=$MPICC
    #ELPA
    if [ $ELPA_MODULES ]
    then
      CLI_ELPA_OPENMP=1
      FLEUR_LIBDIR="$FLEUR_LIBDIR $ELPA_LIB"
      FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $ELPA_MODULES"
    fi
    #hdf5
    if [ $HDF5_ROOT ]
    then 
       FLEUR_LIBDIR="$FLEUR_LIBDIR $HDF5_ROOT/lib"
       FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $HDF5_ROOT/include"
    fi
 
    FLEUR_LIBRARIES="-lelpa_openmp;-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64"
    
