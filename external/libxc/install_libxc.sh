libxc_version=4.2.1
if [ ! -d libxc-${libxc_version} ]
then
    #Get the file with the code
    curl --connect-timeout 10 -LO "https://github.com/MRedies/libxc-mirror/raw/master/libxc-${libxc_version}.tar.gz"
    
    if [ ! -f libxc-${libxc_version}.tar.gz ]; then
        echo "No file found try source"
        curl --connect-timeout 10 -LO "http://www.tddft.org/programs/octopus/download/libxc/${libxc_version}/libxc-${libxc_version}.tar.gz"        
    fi

    tar xzf libxc-${libxc_version}.tar.gz
    cd libxc-${libxc_version}
    #Compile&test (This will take a while)
    ./configure --prefix=$PWD/INSTALL_DIR
    make -j
    make install
else
    cd libxc-${libxc_version}
fi
#Store the installation location
FLEUR_LIBDIR="$PWD/INSTALL_DIR/lib $FLEUR_LIBDIR"
FLEUR_INCLUDEDIR="$PWD/INSTALL_DIR/include $FLEUR_INCLUDEDIR"
if [ $FLEUR_LIBRARIES ]
then
   FLEUR_LIBRARIES="$FLEUR_LIBRARIES;-lxcf03;-lxc"
else
   FLEUR_LIBRARIES="-lxcf03;-lxc"
fi
