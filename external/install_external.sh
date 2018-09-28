#Helper to install external library
#1. creates a external directory in the build location
#2. copies the external/XXX directory for library XXX from source to build
#3. starts the install_XXX.sh script in buildname/external to do the work ....

function compile_external(){
    echo "Installing $library"  
    if [ ! -r external ]
    then
	mkdir external
    fi
    
    cp -r $DIR/external/$library external

    here=$PWD
    
    cd external/$library
    
    . install_${library}.sh

    cd $here
    echo "$library done"
}
