#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

{
function clone_from_git(){
     rm -rf $DIR/configure.sh
    cd $DIR
    git init .
    git remote add -f origin  https://iffgit.fz-juelich.de/fleur/fleur.git
    git checkout release
    echo "Current FLEUR source code obtained from iffgit.fz-juelich.de"
    echo "Restart the configure.sh script"
    exit
    }

function check_git(){
    git_found=`which git`
    ping -q -c1 -W1 iffgit.fz-juelich.de
    iffgit_reachable=$?
    if [ $git_found ] && [ $iffgit_reachable == 0 ] ; then
       return 0;
    fi
    return 1;
}

function update_git(){
    if [ -r $DIR/.git ]
    then
	cd $DIR ; git pull ; cd -
    else
	echo "WARNING........."
	echo "You asked to pull the current version from IFFGIT"
	echo "So far your version is not controlled by git."
	echo "If you modified any source code in your directory it will be deleted"
	echo "Interrupt the script now to abort or press ENTER to continue"
	echo ""
	read y
	clone_from_git
    fi
    }

if [ ! -r $DIR/cmake ]
then
    if check_git
    then
	clone_from_git
    else
	echo "ERROR: No source present"
	echo "You have no 'git' executable and/or the iffgit server is not reachable"
	echo "Please download the complete source manually to your machine"
	exit
    fi
fi

#These contain functions to be used later on...
. $DIR/cmake/machines.sh
. $DIR/external/install_external.sh

#variables to store arguments

machine=""
label=""
backup=0
gitupdate=0
debug=0
error=""
external_lib=""


echo "------------ Welcome to the FLEUR configuration script -------------"

. $DIR/cmake/process_arguments.sh


if [ $gitupdate -gt 0 ]
then
    update_git
fi


#Check if we are using the git version and update if pull was used as an argument
if test -d $DIR/.git
then
    #Check if hook is installed and install it if needed
    if test -h $DIR/.git/hooks/pre-commit
    then
        echo "Git version found"
    else
        mkdir -p $DIR/.git/hooks
        ln -s $DIR/tests/git-hooks/pre-commit $DIR/.git/hooks
        echo "Git version found, hook installed"
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
   echo "Debug version will be built"
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
        mkdir $buildname
    else
	echo "Overwriting old build"
        cd $buildname
        for file in *
        do
          if [[ "$file" == "external" ]] || [[ "$file" == "Testing" ]]
          then
            echo "Keeping $file directory"
          else
            rm -r $file
          fi
        done
        cd -
    fi
else
   mkdir $buildname
fi
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

if [ "$use_ninja" ]
then
    echo "Using Ninja to build"
    NINJAARG="-G Ninja"
else
    NINJAARG=""
fi

# compile external libraries if needed
for library in ${external_lib}
do
    compile_external
done

. $DIR/cmake/store_environment.sh


${cmake} $CMAKE_OPTIONS $NINJAARG -DCMAKE_BUILD_TYPE=$BUILD $DIR 2>&1 |tee configure.out

cd -

if [ ! -r $buildname/Makefile ]
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
}
exit
