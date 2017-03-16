#!/usr/bin/env bash

#This skript is used by buildbot to set up the environment on different machines
# is loads modules and sets environment variables by sourcing the corresponding 
# files in the machines directory and then calls the configuration script with
# the correct machine name

# the environment variable BUILDSLAVE_MACHINE must be set correctly for this to work

m=$BUILDSLAVE_MACHINE


if [[ $m =~ "JURECA-GPU" ]]
then
    source cmake/machines/JURECA/pgisource.sh
    configure.sh JURECA
    exit
fi

if [[ $m =~ "JURECA-GCC" ]]
then
    source cmake/machines/JURECA/gccsource.sh
    configure.sh AUTO
    exit
fi

if [[ $m =~ "JURECA" ]]
then
    source cmake/machines/JURECA/intelsource.sh
    configure.sh JURECA
    exit
fi

configure.sh $m


