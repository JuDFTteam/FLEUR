
find_compilers(){
compiler_list=$1
ccomp=()
ccomp_d=()
for comp in $compiler_list
do
   if  whereis -b $comp 2>&1 >/dev/null
   then
      path=`whereis 2>/dev/null -b $comp|cut -f2 -d:`
      for p in $path
      do
         if [ ! -d $p ]
         then
            if $p --version 2> /dev/null >/dev/null
            then
               desc=`$p --version 2>/dev/null `
               new=true
               if (( ${#desc} < 1 )); then new=false ; fi
               if [[ $comp == mpi* ]];then desc="MPI:$desc" ;fi
               for i in ${!ccomp_d[@]}
               do 
                  olddesc=${ccomp_d[$i]}
                  if [ "$desc" == "$olddesc" ] ;then 
                     new=false
                  fi
               done
               if $new; then
                  ccomp=(${ccomp[@]} $p)
                  IFS=""
                  ccomp_d=(${ccomp_d[@]} "${desc}")
                  IFS=" "
               fi   
            fi
         fi
      done
   fi
done
}

select_compiler(){
echo $1
for i in ${!ccomp[@]}
do
echo -e "${GREEN}$i) ${ccomp[$i]}${NC}"
echo "-------------------------------------------"
echo ${ccomp_d[$i]}
done
read number
if [ $number == '' ];then echo "Invalid choice" ; exit ;fi
if (( $number < 0)) ;then echo "Invalid number" ; exit ;fi
if (( $number >= ${#ccomp[@]} ));then echo "Invalid number" ; exit ;fi
}

interactive_compiler_selection(){
echo "Searching Compilers. Please wait...."
find_compilers "mpif90 mpiifx mpiifort ifx ifort nvfortran gfortran"
select_compiler "Choose the Fortran-Compiler" 
export FC=${ccomp[$number]}

echo "Searching Compilers. Please wait...."
find_compilers "mpicc mpiicc mpiicx icc icx nvc gcc"
select_compiler "Choose the C-Compiler" 
export CC=${ccomp[$number]}

echo "Searching Compilers. Please wait...."
find_compilers "mpic++ mpicxx mpigxx mpiicpc mpiicpx icc icpc icpx g++ nvc++"
select_compiler "Choose the C++-Compiler" 
export CXX=${ccomp[$number]}
}

mpi_wrapper(){
   export FC=${FC:=mpif90}
   export CC=${CC:=mpicc}
   export CXX=${CXX:=mpicxx}
}

intel_old(){
   export FC=${FC:=ifort}
   export CC=${CC:=icc}
   export CXX=${CXX:=icpc}
}

intel(){
   export FC=${FC:=ifx}
   export CC=${CC:=icx}
   export CXX=${CXX:=icpx}
}

nvidia(){
   export FC=${FC:=nvfortran}
   export CC=${CC:=nvc}
   export CXX=${CXX:=nvc++}
}

gfortran(){
   export FC=${FC:=gfortran}
   export CC=${CC:=gcc}
   export CXX=${CXX:=g++}
}

configure_compiler(){
   if [[ $compiler == "nvidia" ]] ; then nvidia ; fi
   if [[ $compiler == "intel" ]] ; then intel ; fi
   if [[ $compiler == "intel_old" ]] ; then intel_old ; fi
   if [[ $compiler == "gfortran" ]] ; then gfortran ; fi
   if [[ $compiler == "mpi" ]] ; then mpi_wrapper ; fi   
   if [[ $compiler == "interactive" ]] ; then interactive_compiler_selection ; fi   

   if [ -z ${FC+x} ] || [ -z ${CC+x} ] || [ -z ${CXX+x} ] ; then
      echo -e "${RED}The compilers should be specified using the FC,CC,CXX environment variables.${NC}"
      echo "Not all of these have been specified either before the call to the configure script"
      echo "or by giving a specific compiler toolchain (see -c option) to the script."
      
      if [[ $compiler == "auto" ]] ; then
         echo "You have choosen the AUTO mode in which cmake will try to determine the compilers."
         echo "GOOD LUCK!"
      else
         echo -e "${RED}You can now use an interactive compiler selection.${NC}"
         echo "Press ENTER to continue or abort the script now"
         read
         interactive_compiler_selection
      fi
   fi

   echo "FC: ${FC}"
   echo "CC: ${CC}"
   echo "CXX:${CXX}"

 }


