version=3.11.4-`uname`-x86_64

#Download if not present...
if [ ! -r cmake-${version} ]
then
    #Get the file with the code
    curl -LO "https://cmake.org/files/v3.11/cmake-${version}.tar.gz"
    tar xzf cmake-${hdf5_version}.tar.gz
fi

#add to path
export PATH=$PWD/cmake-${version}/bin:$PATH
