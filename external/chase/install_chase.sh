#get ChASE does not work yet as the git-version of ChASE is not compatible
#we have a copy of the ChASE directory instead
cd src
if [ ! -r ChASE ]
then 
     git clone https://github.com/ChASE-library/ChASE
#     mv ChASE_archived ChASE
fi
cd ..

#build libchase_fleur

if [ ! -r build ] 
then
  mkdir build
  cd build
  cmake ../src
  make install
  cd ..
fi

#Store path to library for FLEUR
#export FLEUR_LIBDIR="$FLEUR_LIBDIR $PWD/build/FLEUR"
export FLEUR_LIBDIR="$FLEUR_LIBDIR $PWD/build/install/lib64"

