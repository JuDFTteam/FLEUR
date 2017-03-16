#!/usr/bin/env bash

#This skript is used by buildbot to set up the environment on different machines
# is loads modules and sets environment variables by sourcing the corresponding 
# files in the machines directory and then calls the configuration script with
# the correct machine name

# the environment variable BUILDSLAVE_MACHINE must be set correctly for this to work

m=$BUILDSLAVE_MACHINE
export FLEUR_CONFIG_MACHINE=$m

if [[ $m =~ "JURECA-GPU" ]]
then
    source cmake/machines/JURECA/pgisource.sh
    export FLEUR_CONFIG_MACHINE=JURECA
elif [[ $m =~ "JURECA-GCC" ]]
then
    source cmake/machines/JURECA/gccsource.sh
    export FLEUR_CONFIG_MACHINE=AUTO
elif [[ $m =~ "JURECA" ]]
then
    source cmake/machines/JURECA/intelsource.sh
    export FLEUR_CONFIG_MACHINE=JURECA
fi

$*


