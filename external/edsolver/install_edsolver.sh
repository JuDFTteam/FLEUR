

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
   edsolver_version="2019-Jun-10"

   if [ ! -r libEDsolver.${edsolver_version} ]
   then
      if [ -r libEDsolver.${edsolver_version}.tar.gz ]
      then
         tar xzf libEDsolver.${edsolver_version}.tar.gz
         cp default.mk libEDsolver.${edsolver_version}/make/
         cd libEDsolver.${edsolver_version}.tar.gz
         make lib
      else
         echo "libEDsolver not present"
      fi
   else
      cp default.mk libEDsolver.${edsolver_version}/make/
      cd libEDsolver.${edsolver_version}
      if [ ! -r "libEDsolver.a" ]
      then
         make lib
      fi
   fi

   #Store the installation location
   FLEUR_LIBDIR="$PWD/ $FLEUR_LIBDIR"
   FLEUR_INCLUDEDIR="$PWD/ $FLEUR_INCLUDEDIR"
fi