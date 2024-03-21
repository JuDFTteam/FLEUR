#get ChASE does not work yet as the git-version of ChASE is not compatible
#we have a copy of the ChASE directory instead
cd src
if [ ! -r ChASE ]
then 
     git clone https://github.com/SimLabQuantumMaterials/ChASE.git
fi
cd ..

#build libchase_fleur

if [ ! -r build ] 
then
  mkdir build
  cd build
  cmake ../src
  make
  cd ..
fi

#Store path to library for FLEUR
export FLEUR_LIBDIR="$FLEUR_LIBDIR $PWD/build/FLEUR"

