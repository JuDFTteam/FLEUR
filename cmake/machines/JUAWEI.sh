  if  module list 2>&1 | grep -q mpich
  then
     export FC=mpif90
     export CC=mpicc
     export CXX=mpicxx
  else 
     export FC=gfortran
     export CC=gcc
     export CXX=g++
  fi
