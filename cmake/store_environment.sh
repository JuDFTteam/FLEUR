# script to store environment settings for later use in cmake

echo "set(compile_user $USER)" >config.cmake

#Some freqeuntly used Environment variables
if [ -n "$HDF5_ROOT" ] ; then FLEUR_LIBDIR="$FLEUR_LIBDIR $HDF5_ROOT/lib" ; FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $HDF5_ROOT/include" ; fi
if [ -n "$HDF5_DIR"} ] ; then FLEUR_LIBDIR="$FLEUR_LIBDIR $HDF5_DIR/lib" ; FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $HDF5_DIR/include" ; fi
if [ -n "$HDF5_LIB" ] ; then FLEUR_LIBDIR="$FLEUR_LIBDIR $HDF5_LIB" ; fi
if [ -n "$HDF5_INCLUDE" ] ; then FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $HDF5_INCLUDE" ; fi
if [ -n "$HDF5_MODULES" ] ; then FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $HDF5_MODULES" ; fi



#Set options for linker
#1. if environment variable FLEUR_LIBRARIES is present use it 
#2. if CLI_LIBRARIES is present, use it
#3. if FLEUR_LIBDIR is present, add these directories with -L option

if [ "$FLEUR_LIBRARIES" ]
then
    cmake_lib="$FLEUR_LIBRARIES"
fi
if [ "$CLI_LIBRARIES" ]
then
    if [ "cmake_lib" ]
    then
	cmake_lib="$CLI_LIBRARIES;$cmake_lib"
    else
	cmake_lib="$CLI_LIBRARIES"
    fi
fi
#check the FLEUR_LIBDIR variable

for lib in $FLEUR_LIBDIR $CLI_LIBDIR 
do
    if [ "cmake_lib" ]
    then
	cmake_lib="-L$lib;$cmake_lib"
    else
	cmake_lib="-L$lib"
    fi
done

echo "set(FLEUR_LIBRARIES $cmake_lib)" >>config.cmake

#Set compiler flags
#1. If CMAKE_Fortran_FLAGS is given use that
#2. If CLI_FLAGS is given use that
#3. If FLEUR_INCLUDEDIR/CLI_INCLUDEDIR is given, add -I options for these

if [ "$CMAKE_Fortran_FLAGS" ]
then
    cmake_flags="$CMAKE_Fortran_FLAGS"
fi
if [ "$CLI_FLAGS" ]
then
    cmake_flags="$CMAKE_Fortran_FLAGS $cmake_flags"
fi
for lib in $FLEUR_INCLUDEDIR $CLI_INCLUDEDIR
do
    cmake_flags="-I$lib $cmake_flags"
done
echo "set(CMAKE_Fortran_FLAGS \"$cmake_flags\")" >>config.cmake

#Set options to turn on/off specific features


if [ "$FLEUR_USE_HDF5" ] || [ "$FLEUR_USE_SERIAL" ] || [ "$FLEUR_USE_MPI" ] || [ "$FLEUR_USE_WANNIER" ] || [ "$FLEUR_USE_MAGMA" ] || [ "$FLEUR_USE_GPU" ]
then
    echo "WARNING"
    echo "The FLEUR_USE_XXX environment variables are no longer supported, use the command line options instead"
fi    

if [ "$CLI_USE_HDF5" ]
then
    echo "set(CLI_FLEUR_USE_HDF5 $CLI_USE_HDF5)"  >>config.cmake
fi

if [ "$CLI_USE_MPI" ]
then
    echo "set(CLI_FLEUR_USE_MPI $CLI_USE_MPI)"  >>config.cmake
fi

if [ "$CLI_USE_WANNIER" ]
then
    echo "set(CLI_FLEUR_USE_WANNIER $CLI_USE_WANNIER)"  >>config.cmake
fi

if [ "$CLI_USE_CHASE" ]
then
    echo "set(CLI_FLEUR_USE_CHASE $CLI_USE_CHASE)"  >>config.cmake
fi

if [ "$CLI_USE_MAGMA" ]
then
    echo "set(CLI_FLEUR_USE_MAGMA $CLI_USE_MAGMA)"  >>config.cmake
fi

if [ "$CLI_USE_GPU" ]
then
    echo "set(CLI_FLEUR_USE_GPU $CLI_USE_GPU)"  >>config.cmake
fi

if [ "$CLI_USE_SERIAL" ]
then
    echo "set(FLEUR_USE_SERIAL $CLI_USE_SERIAL)"  >>config.cmake
else
    echo "set(FLEUR_USE_SERIAL TRUE)"  >>config.cmake
fi

if [ "$CLI_ELPA_OPENMP" ]
then
    echo "set(CLI_ELPA_OPENMP 1)"  >>config.cmake
fi

if [ "$CLI_WARN_ONLY" ]
then
    echo "set(CLI_WARN_ONLY 1)"  >>config.cmake
fi



