  if  module list 2>&1 | grep -q mpt
  then
     export FC=mpif90
     export CC=mpicc
  fi
