  if  module list 2>&1 | grep -q mpich
  then
     export FC=mpif90
     export CC=mpicc
  else 
     export FC=gfortran
     export CC=gcc
  fi
