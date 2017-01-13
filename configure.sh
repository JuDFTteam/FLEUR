#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "------------ Welcome to the FLEUR configuration script -------------"
. $DIR/cmake/machines.sh
#check if -h or  --help was given as argument
if [ "$1" = "" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
   echo "USAGE: configure.sh MACHINE [debug]"
   echo "
  To help the script finding a proper configuration you should
  provide the name of a specific machine to compile on.
  Currently known machine configurations are:
  " 
   echo "   $known_machines"
   echo " 
  If you do not find a proper choice you might try
        'configuration.sh AUTO'

  You might also want to add your configuration to the file 
  cmake/machines.sh in this case :-)
  
  In addition you can modify some environment variables:
        FLEUR_NO_SERIAL     -- if defined no serial executables will be build
        FC                  -- name of Fortran compiler
        CC                  -- name of C compiler
        FLEUR_LIBRARIES     -- list of linker arguments i.e. '-L/lib;-lbla'
        CMAKE_Fortran_FLAGS -- list of compiler options i.e. '-r8'"
echo "
   By specifying 'debug' in addition to your machine configuration you will build a debugging version"
fi
#Check if we are using the git version and ask if we want to update
if test -d $DIR/.git
then
   #We are using the git version so ask the user (for 10sec)
   echo "Shall we try to update to the newest git version? (y/n)"
   read -n 1 -t 10 x
   if test "$x" == "y"
   then
       cd $DIR 
       git pull
       cd -
   fi
fi


#Now check the machine and set some defaults 
machine=$1
configure_machine

#include a configfile if present
if test -r config.sh
then
 . config.sh
fi

#check if there is a build directory
if test -d build
then
    echo "OLD build directory found, saved in build.$$"
    mv build build.$$
fi
mkdir build
cd build
#run cmake
if test "debug" == "$2"
then
   echo "Debug version will be build"
   BUILD=Debug
else
   BUILD=Release
fi
cmake -DCMAKE_BUILD_TYPE=$BUILD $DIR


echo "Configuration finished"
echo "If no errors occured you should change into directory 'build' "
echo "run 'make' or 'make -j4'"
