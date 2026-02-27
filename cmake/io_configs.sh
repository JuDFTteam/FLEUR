store_config(){
    filename=".$label.sh"
    echo "#!/bin/sh"> $filename
    #CHECK for modules
    if `which module  >/dev/null 2>&1` 
    then
        module save $label
        echo "module restore $label" >>$filename
    fi
    #Compilers
    if [[ ! -z ${CF} ]] ; then echo "export CF=$CF" >> $filename ;fi
    if [[ ! -z ${CC} ]] ; then echo "export CC=$CC" >> $filename ;fi
    if [[ ! -z ${CXX} ]] ; then echo "export CXX=$CXX" >> $filename ;fi
    #FLAGS
    if [[ ! -z ${CFLAGS} ]] ; then echo "export CFLAGS=$CFLAGS" >> $filename ;fi
    if [[ ! -z ${FFLAGS} ]] ; then echo "export FFLAGS=$FFLAGS" >> $filename ;fi
    if [[ ! -z ${CXXFLAGS} ]] ; then echo "export CXXFLAGS=$CXXFLAGS" >> $filename ;fi
    if [[ ! -z ${CPPFLAGS} ]] ; then echo "export CPPFLAGS=$CPPFLAGS" >> $filename ;fi
    if [[ ! -z ${LDFLAGS} ]] ; then echo "export CPPFLAGS=$LDFLAGS" >> $filename ;fi
    #Put call into config-file    
    echo "$command" >>$filename
    chmod +x $filename
}