#This file contains the defaults for compiling FLEUR on known machines

#please add further definitions here and also add code below!
read -r -d '' known_machines << EOM
   IFF         -- cluster @ PGI of FZJ
   JURECA      -- @JSC
   JUQUEEN     -- @JSC
   CLAIX       -- @RWTH
   MARCONI     -- @CINECA 
EOM


function configure_machine(){
    
    if [ "$machine" = "JURECA" ] 
    then
	echo "JURECA configuration used"
	if module list 2>&1 |grep -q -i intel
	then
	    echo "Intel toolchain used"
	    if module list 2>&1| grep -q Python &&
               module list 2>&1| grep -q CMake &&
	       module list 2>&1| grep -q ELPA
               #module list 2>&1| grep -q HDF5 
	    then
		echo "All required modules load loaded"
	    else
		echo "You have to load the required modules"
		echo "See and/or source $DIR/cmake/machines/JURECA/intelsource.sh"
		exit
	    fi
	    cp $DIR/cmake/machines/JURECA/JURECA.cmake config.cmake
	elif module list 2>&1 |grep -q PGI
	then
	    echo "PGI toolchain used"
	    if module list 2>&1| grep -q MVAPICH2 &&
               module list 2>&1| grep -q CMake &&
	       [ -n "${MAGMA_ROOT}" ] &&
	       [ -n "${XML2_ROOT}" ] 
            then
		echo "All required modules loaded, all variables set"
	    else
		echo "Not all modules are loaded"
		echo "See and/or source $DIR/cmake/machines/JURECA/pgisource.sh"
		echo "And set the variables XML2_ROOT and MAGMA_ROOT"
		exit
	    fi
	    cp $DIR/cmake/machines/JURECA/JURECAGPU.cmake config.cmake
	else
	    echo "You have to load the correct modules for compiling"
	    echo " Look for files to source in $DIR/cmake/machines/JURECA"
	    exit
	fi
    # JUQUEEN
    elif [ "$machine" = "JUQUEEN" ]
    then 
   	echo "JUQUEEN configuration used"
	if module list 2>&1| grep -q hdf5 &&
           module list 2>&1| grep -q scalapack
        then
           echo "All required modules load loaded"
	else
	   echo "You have to load the required modules"
	   echo "module load hdf5/1.8.15_BGQ scalapack/2.0.2_elpa_simd"
	   exit
	fi
        cp $DIR/cmake/machines/JUQUEEN/JUQUEEN.cmake config.cmake
	
    #IFF linux cluster
    elif [ "$machine" = "IFF" ]
    then
	echo "IFF cluster configuration used"
	cp $DIR/cmake/machines/IFF.cmake config.cmake

    #RWTH cluster
    elif [ "$machine" = "CLAIX" ]
    then
	echo "CLAIX@RWTH configuration used"
	if ! module list 2>&1| grep -q intelmpi
	then
	    echo "Please use intelmpi, e.g. do a module switch openmpi intelmpi"
	    exit
	fi
	cp $DIR/cmake/machines/CLAIX.cmake config.cmake
module load LIBRARIES
    elif [ "$machine" = "MARCONI" ]
    then
	if ! module list 2>&1| grep -q " intel\/" || ! module list 2>&1| grep -q " intelmpi" ||! module list 2>&1| grep -q cmake
	then
	    echo "Load the modules needed to compile: intel,intelmpi,cmake"
	    exit
	fi
	cp $DIR/cmake/machines/MARCONI.cmake config.cmake
    elif [ "$machine" = "AUTO" ] 
    then
	echo "No machine specific settings used"
	echo "GOOD LUCK!"
    else
	echo "No valid machine configuration specified"
	exit
    fi
}

