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
    if ! echo "$module_list" |grep ELPA | grep hybrid
    then
         error=1
    fi     
    if ! echo "$module_list" |grep HDF5
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
    #hdf5
    FLEUR_LIBDIR="$FLEUR_LIBDIR $HDF5_DIR/lib"
    FLEUR_INCLUDEDIR="$FLEUR_INCLUDEDIR $HDF5_DIR/include"
 
 
    FLEUR_LIBRARIES="-lelpa_openmp;-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64"
    
elif module list 2>&1 |grep -q PGI
then
    echo "PGI toolchain used"
    FC=mpif90
    CC=mpicc
    FLEUR_INCLUDEDIR="$XML2_ROOT/include"
    FLEUR_LIBRARIES="-lstdc++;-L$MKLROOT/lib/intel64;-lmkl_scalapack_lp64;-lmkl_intel_lp64;-lmkl_pgi_thread;-lmkl_core;-lmkl_blacs_intelmpi_lp64"
else
    echo "You need to load the modules for the compiler"
    echo "e.g. module load intel-para"
    exit
fi
