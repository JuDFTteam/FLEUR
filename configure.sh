#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#process arguments
help=0
all_tests=0
machine=""
label=""
backup=0
gitupdate=0
debug=0
error=""
while [ $# -gt 0 ]
do
    case "$1" in
        -h) help=1;;
	-help) help=1;;
        -b) backup=1;;
        -backup) backup=1;;
        -g) gitupdate=1;;
	-gitupdate) gitupdate=1;;
	-t) all_tests=1;;
	-all_tests) all_tests=1;;
	-l) shift;label=$1;;
	-m) shift;machine=$1;;
	-*) error="Unkown argument";;
	*)  break;;	# terminate while loop
    esac
    shift
done
if [ $# -gt 0 ]
then
    if [ "$machine" = "" ]
    then
	machine=$1;shift
    else
	error="You specified the -m switch and gave an additional MACHINE argument" 
    fi
fi
if [ $# -gt 0 ]
then
    if [ "$label" = "" ]
    then
	label=$1;shift
    else
	error="You specified the -l switch and gave an additional LABEL argument"
    fi
fi
if [ $# -gt 0 ]
then
    error="Extra unkown arguments"
fi


echo "------------ Welcome to the FLEUR configuration script -------------"

if [ "$error" != "" ]
then
    echo $error
    echo "ERROR in calling configure. Please use -h for help on usage"
    exit 1
fi

. $DIR/cmake/machines.sh

#check if -h or  -help was given as argument
if [ $help -gt 0 ] 
then
   echo "USAGE: configure.sh [options] [MACHINE] [label]"
   echo "
  The following command line options are supported:
  -h   : print this help-page
  -m # : specify the machine to build on (see below for possible options)
         This can also be specified without -m as a first argument
  -l # : label for the build. It will be attached to the name of the build directory.
         This option can also be specified as second argument without -l
  -d   : build a debugging version of FLEUR (adds .debug to label)
  -g   : do a git pull first if this is a git version
  -t   : generate all tests including those that run longer
  -b   : backup an old build directory if present
  
  To help the script finding a proper configuration you should
  provide the name of a specific machine to compile on.
  Currently known machine configurations are:
  " 
   echo "   $known_machines"
   echo " 
  If you do not specify the machine the AUTO option will be used in which some
  defaults are tested. It might work if your machine is not one of those known.

  You might also want to add your configuration to the file 
  cmake/machines.sh in this case :-)
  
  In addition you can modify some environment variables:
        FC                  -- name of Fortran compiler
        CC                  -- name of C compiler
        FLEUR_LIBRARIES     -- list of linker arguments i.e. '-L/lib;-lbla'
        CMAKE_Fortran_FLAGS -- list of compiler options i.e. '-r8'"
echo "
   By specifying a label you can have different build directories.
   The label will be added to the 'build' directory name."
  exit 1
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
   if [ $gitupdate -gt 0 ] 
   then
       cd $DIR 
       git pull
       cd -
   fi
fi


#include a configfile if present
if test -r config.sh
then
 . config.sh
fi

#Name of the build directory
if [ -n "$label" ]
then
    buildname="build.$label"
else
    buildname="build"
fi

if [ $debug -gt 0 ]
then
   echo "Debug version will be build"
   BUILD=Debug
   $buildname="$buildname.debug"
else
   BUILD=Release
fi

#check if there is a build directory
if test -d $buildname
then
    if [ $backup -gt 0 ]
    then
	echo "OLD build directory found, saved in build.$$"
	mv $buildname $buildname.$$
    else
	echo "Overwriting old build"
	rm -r $buildname
    fi  
fi
mkdir $buildname
cd $buildname

#Now check the machine and set some defaults 
if [[ $machine =~ FLEUR_CONFIG_MACHINE ]]
then
    machine=$FLEUR_CONFIG_MACHINE
fi
if [ "$machine" = "" ]
then
    machine=AUTO
fi
echo "Machine config:$machine"
configure_machine

#run cmake

cmake -Dall_tests=$all_tests -DCMAKE_BUILD_TYPE=$BUILD $DIR

if [ -r $buildname/Makefile ]
then
    echo "Your configuration failed"
    echo "Perhaps you have to specify compiler options or give a machine dependent configuration."
    echo "You might want to call the configure.sh script with -h"
else
    echo "Configuration finished"
    echo "You should change into directory $buildname "
    echo "run 'make' or 'make -j'"
fi
