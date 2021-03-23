# Configuration for JURECA-DC, GCC compiler

#determine XML2 module
xml2=$(module --terse --show_hidden spider libxml2 2>&1 | grep -v '^$' | tail -n 1)

module_list=$(module --terse list 2>&1)

if echo "$module_list" |grep -q -i GCC
then
    echo "GCC toolchain used."

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
	exit
    fi

    export FC=${FC:=mpif90}
    export CC=${CC:=mpicc}
    export CXX=${CXX:=mpicxx}
    #ELPA
    CLI_ELPA_OPENMP=1
    FLEUR_LIBDIR="$FLEUR_LIBDIR $ELPA_LIB"
    FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $ELPA_MODULES_OPENMP"
else
    echo "You need to load the modules for the compiler"
    exit
fi
