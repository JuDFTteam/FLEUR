elpa_version=2017.11.001
if [ ! -r elpa-${elpa_version} ]
then
    #Get the file with the code
    curl -LO "http://elpa.mpcdf.mpg.de/html/Releases/${elpa_version}/elpa-${elpa_version}.tar.gz"
    tar xzf elpa-${elpa_version}.tar.gz
    cd elpa-${elpa_version}
    #configure
    ./configure --disable-doxygen-doc --enable-shared=no --disable-mpi-module --enable-openmp --prefix=$PWD/INSTALL_DIR
    #Compile&test (This will take a while)
    make
    #Do make install 
    make install
    cd ..
fi
#Store the installation location
FLEUR_LIBDIR="$PWD/INSTALL_DIR/lib $FLEUR_LIBDIR"
FLEUR_INCLUDEDIR="$PWD/INSTALL_DIR/include $FLEUR_INCLUDEDIR"
