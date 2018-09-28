#Configuration for JURECA@FZJ; Use intel-toolchain, e.g. do a 'module load intel-para'
#determine XML2 module
xml2=`module --show_hidden spider libxml2 2>&1 |grep libxml2/|grep -v module |tail -1`
#determine ELPA module
elpa=`module spider ELPA 2>&1 |grep hybrid`


if module list 2>&1 |grep -q -i intel
then
    echo "Intel toolchain used"

    #check if all modules are loaded
    module_list=`module list 2>&1`
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
        #determine XML2 module
        xml2=`module --show_hidden spider libxml2 2>&1 |grep libxml2/|grep -v module |tail -1`
        #determine ELPA module
        elpa=`module spider ELPA 2>&1 |grep hybrid`
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
    echo "e.g. module load intel-para CMake HDF5 $xml2 $elpa "
    exit
fi
