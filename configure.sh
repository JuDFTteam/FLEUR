#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#These contain functions to be used later on...
. $DIR/cmake/machines.sh
.  $DIR/external/install_external.sh

#variables to store arguments

all_tests=0
machine=""
label=""
backup=0
gitupdate=0
debug=0
error=""
external_lib=""


echo "------------ Welcome to the FLEUR configuration script -------------"

. $DIR/cmake/process_arguments.sh



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
   buildname="$buildname.debug"
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
if [ "$cmake" ]
then
    echo "Using provided cmake:$cmake"
else
    cmake="cmake"
fi

# compile external libraries if needed
for library in ${external_lib}
do
    compile_external
done

. $DIR/cmake/store_environment.sh


${cmake} $CMAKE_OPTIONS -Dall_tests=$all_tests -DCMAKE_BUILD_TYPE=$BUILD $DIR 2>&1 |tee configure.out

if [ -r $buildname/Makefile ]
then
    echo "Your configuration failed"
    echo "Perhaps you have to specify compiler options or give a machine dependent configuration."
    echo "You might want to call the configure.sh script with -h"
else
    echo "Configuration finished"
    if [ "$make_directly" ]
    then
	cd $buildname
	make -j
    else
	echo "You should change into directory $buildname "
	echo "run 'make' or 'make -j'"
    fi
fi
