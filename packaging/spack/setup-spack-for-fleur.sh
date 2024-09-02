find_spack(){
    SPACK=${SPACK:-`which spack`}

    if ! which "$SPACK" >/dev/null
    then
        echo " ========= ERROR ======================  "
        echo -e "SPACK not found"
        echo "To use this function, please install spack as described on"
        echo "https://spack.io/about"
        echo ""
        echo "After activating the spack environement using e.g.:"
        echo "source $SPACKDIR/share/spack/setup-env.sh"
        echo "try again."
        exit
    fi
}


create_spack_env(){
    echo 
    sed -e "s@..FLEUR_SOURCE.@${DIR}/@" $env_file >tmp.yaml
    $SPACK env create $newname tmp.yaml
    rm tmp.yaml
}

get_environments(){
    for e in ${FLEURDIR}/environments/*.yaml
    do
        spack_envs+=($e)
    done    
}

select_env(){
    get_environments
    for i in "${!spack_envs[@]}"; do 
        echo -e "${GREEN}$i\t${spack_envs[$i]}${NC}"
        d=`grep "##docu" ${spack_envs[$i]}|sed -e "s/##docu/      /g"`
        echo -e "$d"
    done
    read -p "Select number of predefined settings to use:" REPLY
    env_file=${spack_envs[$REPLY]}
    read -p "Please enter name of new environment:" newname
    create_spack_env
}

setup_spack(){
    FLEURDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
    find_spack
    select_env
}