if  module list 2>&1 | grep -q mpt
then
   export FC=mpif90
   export CC=mpicc
   export CXX=mpicxx
   FLEUR_LIBRARIES="-lxml2;-lmkl_scalapack_lp64;-lmkl_blacs_sgimpt_lp64"
   #FLEUR_LIBRARIES="-lxml2"
elif  module list 2>&1 | grep -q openmpi
then
   export FC=mpif90
   export CC=mpicc
   FLEUR_LIBRARIES="-lxml2;-lmkl_scalapack_lp64;-lmkl_blacs_openmpi_lp64"
elif  module list 2>&1 | grep -q impi
then
   #export FC=mpif90
   #export CC=mpicc
   export FC=mpiifort
   export CC=mpiicc
   export CXX=mpiicpc
   FLEUR_LIBRARIES="-lxml2;-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64"
fi

#if  module list 2>&1 | grep -q elpa
#then
#   #ELPA
#   CLI_ELPA_OPENMP=1
#   FLEUR_LIBDIR="$FLEUR_LIBDIR "
#   FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $ELPA_MODULES_OPENMP"
#fi

echo "Wannier is switched off-manually"
CLI_USE_WANNIER=FALSE
echo "ChASE is switched off-manually"
CLI_USE_CHASE=FALSE

#ELPA
if [ $ELPA_MODULES ]
then
  CLI_ELPA_OPENMP=1
  FLEUR_LIBDIR="$FLEUR_LIBDIR $ELPA_LIB"
  FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $ELPA_MODULES"
fi


