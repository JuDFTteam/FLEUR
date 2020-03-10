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
echo "------------ Welcome to the FLEUR configuration script -------------"
echo "  This version of the script will only attempt to download the code "
echo "  Please put it in the directory which should contain the source    "
echo "--------------------------------------------------------------------"
if [ $DIR != $PWD ]
then
    echo "WARNING: code will not be downloaded into your current working dir: $PWD"
    echo "Instead it will be put to the directory of the script:$DIR"
    echo "Press ENTER to continue"
    read y
fi
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
}

exit
