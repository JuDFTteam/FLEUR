# Configuration for JURECA/JUWELS@FZJ; Use intel-toolchain, e.g. do a 'ml intel-para'

#determine XML2 module
xml2=$(module --terse --show_hidden spider libxml2 2>&1 | grep -v '^$' | tail -n 1)
#determine ELPA module
elpa=`module spider ELPA 2>&1 |grep hybrid`

module_list=$(module --terse list 2>&1)

if echo "$module_list" |grep -q -i intel
then
    echo "Intel toolchain used"

    #check if all modules are loaded
    error=0
    if ! echo "$module_list" |grep CMake
    then
         error=1
    fi     
    if ! echo "$module_list" |grep xml2
    then
         error=1
    fi     
    if [[ "$error" -gt 0 ]]
    then
	echo "Not all required modules are loaded"
	echo "You should do:"
    	echo
	echo "ml intel-para CMake HDF5 $xml2 $elpa"
	echo
	echo "It might be a good idea to include that in your profile"
	exit
    fi

    export FC=${FC:=mpif90}
    export CC=${CC:=mpicc}
    #ELPA
    CLI_ELPA_OPENMP=1
    FLEUR_LIBDIR="$FLEUR_LIBDIR $ELPA_LIB"
    FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $ELPA_MODULES_OPENMP"
else
    echo "You need to load the modules for the compiler"
    echo "e.g. ml intel-para CMake HDF5 $xml2 $elpa "
    exit
fi
