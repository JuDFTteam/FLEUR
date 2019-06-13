libxml2_version=2.7.2
if [ ! -r libxml2-${libxml2_version} ]
then
    #Get the file with the code
    curl -LO "ftp://xmlsoft.org/libxml2/libxml2-${libxml2_version}.tar.gz"
   
    tar xzf libxml2-${libxml2_version}.tar.gz
    cd libxml2-${libxml2_version}
    #Compile&test (This will take a while)
    ./configure  --disable-shared --without-python
    make 
else
    cd libxml2-${libxml2_version}    
fi
#Store the installation location
FLEUR_LIBDIR="$PWD/.libs $FLEUR_LIBDIR"
FLEUR_INCLUDEDIR="$PWD/include $FLEUR_INCLUDEDIR"
