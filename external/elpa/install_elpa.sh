#elpa_version=2017.11.001
elpa_version=2018.05.001.rc1
if [ ! -r elpa-${elpa_version} ]
then
    #Get the file with the code
    curl -LO "http://elpa.mpcdf.mpg.de/html/Releases/${elpa_version}/elpa-${elpa_version}.tar.gz"
    tar xzf elpa-${elpa_version}.tar.gz
    cd elpa-${elpa_version}
    #configue
    ./configure --enable-c-tests=no --disable-doxygen-doc --enable-shared=no --disable-mpi-module --enable-openmp --prefix=$PWD/INSTALL_DIR \
                --enable-option-checking=fatal SCALAPACK_LDFLAGS="-L$MKLROOT/lib/intel64 -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread "\
                SCALAPACK_FCFLAGS="-I$MKL_HOME/include/intel64/lp64"
    #Compile&test (This will take a while)
    make
    #Do make install 
    make install
    cd ..
fi
#Store the installation location
if [[ -r $PWD/INSTALL_DIR/lib/libelpa.a ]] || [[ -r $PWD/INSTALL_DIR/lib/libelpa_openmp.a ]]
then
      echo "libelpa.a found"
else	
    #Try to copy files manually if make install failed
    mkdir $PWD/INSTALL_DIR
    mkdir $PWD/INSTALL_DIR/lib/
    cp elpa-${elpa_version}/.libs/* $PWD/INSTALL_DIR/lib/
    mkdir $PWD/INSTALL_DIR/include
    cp elpa-${elpa_version}/modules/* elpa-${elpa_version}/private_modules/* $PWD/INSTALL_DIR/include
fi
FLEUR_LIBDIR="$PWD/INSTALL_DIR/lib $FLEUR_LIBDIR"
FLEUR_INCLUDEDIR="$PWD/INSTALL_DIR/include $FLEUR_INCLUDEDIR"
if [ -r $PWD/INSTALL_DIR/lib/libelpa_openmp.a ]
then
CLI_ELPA_OPENMP=1
fi 
