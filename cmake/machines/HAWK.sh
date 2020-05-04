if  module list 2>&1 | grep -q mpt
then
   export FC=mpif90
   export CC=mpicc
fi

FLEUR_LIBRARIES="-lxml2;-lmkl_scalapack_lp64;-lmkl_blacs_sgimpt_lp64"
echo "Wannier is switched off-manually"
CLI_USE_WANNIER=FALSE
echo "ChASE is switched off-manually"
CLI_USE_CHASE=FALSE


