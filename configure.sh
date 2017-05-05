#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "------------ Welcome to the FLEUR configuration script -------------"
. $DIR/cmake/machines.sh
#check if -h or  --help was given as argument
if [ "$1" = "" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
   echo "USAGE: configure.sh MACHINE [label]"
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
        FC                  -- name of Fortran compiler
        CC                  -- name of C compiler
        FLEUR_LIBRARIES     -- list of linker arguments i.e. '-L/lib;-lbla'
        CMAKE_Fortran_FLAGS -- list of compiler options i.e. '-r8'"
echo "
   By specifying a label which contains 'debug' in addition to your 
   machine configuration you will build a debugging version. Otherwise
   the label will be added to the build directory name."
fi

#Check if we are using the git version and update if pull was used as an argument
if test -d $DIR/.git 
    #Check if hook is installed and install it if needed
    if test -h $DIR/.git/hooks/pre-commit
    then
        echo "Git version found"
    else
        ln -s $DIR/tests/git-hooks/pre-commit $DIR/.git/hooks
        echo "Git version found, hook installed"
    fi
then
   if [[ $1 =~ .*pull.* ]] || [[ $2 =~ .*pull.* ]] || [[ $3 =~ .*pull.* ]] 
   then
       cd $DIR 
       git pull
       cd -
       exit
   fi
fi


#include a configfile if present
if test -r config.sh
then
 . config.sh
fi

#Name of the build directory
label=$2
if [ -n "$label" ]
then
    buildname="build.$label"
else
    buildname="build"
fi

#check if there is a build directory
if test -d $buildname
then
    echo "OLD build directory found, saved in build.$$"
    mv $buildname $buildname.$$
fi
mkdir $buildname
cd $buildname

#Now check the machine and set some defaults 
machine=$1
if [[ $machine =~ FLEUR_CONFIG_MACHINE ]]
then
    machine=$FLEUR_CONFIG_MACHINE
fi
echo "Machine config:$machine"
configure_machine

#run cmake
if [[ $buildname =~ .*debug.* ]]
then
   echo "Debug version will be build"
   BUILD=Debug
else
   BUILD=Release
fi
cmake -DCMAKE_BUILD_TYPE=$BUILD $DIR


echo "Configuration finished"
echo "If no errors occured you should change into directory $buildname "
echo "run 'make' or 'make -j4'"
