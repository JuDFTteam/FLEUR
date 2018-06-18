libxc_version=4.2.1
#Get the file with the code
curl -LO "http://www.tddft.org/programs/octopus/download/libxc/${libxc_version}/libxc-${libxc_version}.tar.gz"
tar xzf libxc-${libxc_version}.tar.gz
cd libxc-${libxc_version}
#Compile&test (This will take a while)
./configure --prefix=$PWD/INSTALL_DIR
make
make install

#Store the installation location
FLEUR_LIBDIR="$PWD/INSTALL_DIR/lib $FLEUR_LIBDIR"
FLEUR_INCLUDEDIR="$PWD/INSTALL_DIR/include"
if [ $FLEUR_LIBRARIES ]
then
   FLEUR_LIBRARIES="$FLEUR_LIBRARIES;-lxcf03;-lxc"
else
   FLEUR_LIBRARIES="-lxcf03;-lxc"
fi
