FC=${FC:=mpiifort}
FLEUR_LIBRARIES="-lxml2;-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64;-mt_mpi;${FLEUR_LIBRARIES}"
