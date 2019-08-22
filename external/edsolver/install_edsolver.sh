edsolver_version="2019-Aug-22"

if [ "$machine" = "AUTO" ]
then
   ED_PLATFORM="gfortran"
else
   if [ "$machine" = "INTEL" ] || [ "$machine" = "INTEL_MPI" ]
   then
      ED_PLATFORM="ifort"
   else
      echo "libEDsolver at the moment only configured with AUTO and INTEL/INTEL_MPI"
   fi
fi


#Add the linking arguments for fleur for both the EDsolver and ARPACk
CLI_LIBRARIES="-lEDsolver;-larpack_$ED_PLATFORM $CLI_LIBRARIES"

if [ ! -r  libEDsolver_FLEUR ]
then 
   if [ -r libEDsolver_FLEUR.tar.gz ]
   then
      tar xzf libEDsolver_FLEUR.tar.gz
   else
      echo "No valid version of EDsolver library packaged for FLEUR available"
   fi
fi

if [ -r libEDsolver_FLEUR ]
then
   cd libEDsolver_FLEUR
   #First add ARPACK

   if [ ! -r ARPACK ]
   then
      echo "ARPACK is needed for the EDsolver"
   else
      FLEUR_LIBDIR="$PWD/ARPACK $FLEUR_LIBDIR"
   fi


   #Now look at the EDsolver

   if [ ! -r libEDsolver.${edsolver_version} ]
   then
      if [ -r libEDsolver.${edsolver_version}.tar.gz ]
      then
         tar xzf libEDsolver.${edsolver_version}.tar.gz
         cd libEDsolver.${edsolver_version}.tar.gz
         make lib PLATFORM="$ED_PLATFORM"
      else
         echo "libEDsolver not present"
      fi
   else
      cd libEDsolver.${edsolver_version}
      if [ ! -r "libEDsolver.a" ]
      then
         make lib PLATFORM="$ED_PLATFORM"
      fi
   fi

   #Store the installation location
   FLEUR_LIBDIR="$PWD/ $FLEUR_LIBDIR"
   FLEUR_INCLUDEDIR="$PWD/ $FLEUR_INCLUDEDIR"

   #Make the utils (if the compilation was sucessful)
   if [ -r "libEDsolver.a" ]
   then
      cd utils
      make all PLATFORM="$ED_PLATFORM"
   fi
fi