#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#color for output
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

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
	echo -e "${RED}WARNING.........${NC}"
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
	echo -e "${RED}ERROR: No source present${NC}"
	echo "You have no 'git' executable and/or the iffgit server is not reachable"
	echo "Please download the complete source manually to your machine"
	exit
    fi
fi



#variables to store arguments

compiler=""
label=""
backup=0
gitupdate=0
debug=0
error=""
external_lib=""


echo -e "${RED}------------ Welcome to the FLEUR configuration script -------------${NC}"

. $DIR/cmake/process_arguments.sh


if [ $conf_spack -gt 0 ]
then
     source $DIR/packaging/spack/setup-spack-for-fleur.sh
     setup_spack
     exit
fi


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

#Now check the compiler settings

#These contain functions to be used now ...
. $DIR/cmake/compiler.sh

compiler=${compiler:=none}
configure_compiler

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
    echo -e "${RED}Your configuration failed"
    echo "Perhaps you have to specify compiler options or choose a different compiler."
    echo -e "You might want to call the configure.sh script with -h${NC}"
else
    echo "Configuration finished"
    if [ "$make_directly" ]
    then
	cd $buildname
	make -j
    else
	echo -e "${GREEN}You should change into directory $buildname "
	echo -e "run 'make' or 'make -j'${NC}"
    fi
fi
}
exit
