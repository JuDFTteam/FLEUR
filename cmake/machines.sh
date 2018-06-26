#This file contains the defaults for compiling FLEUR on known machines

cd $DIR/cmake/machines/
for f in *.sh
do
  config=`basename $f .sh`
  desc=`head -1 $f|cut -d# -f 2`
  if [ ! "$desc" == "NOSHOW" ]
  then
  known_machines="$known_machines
  $config
       -- $desc "
  fi
done
cd -

function configure_machine(){
    if [ "$machine" = "AUTO" ] 
    then
     	echo "No machine specific settings used"
	echo "GOOD LUCK!"
    else
	if [ -r $DIR/cmake/machines/$machine.sh ]
	then
	    echo "Using config for: $machine"
	    . $DIR/cmake/machines/$machine.sh
	else
	    echo "No valid machine configuration specified"
	    exit
	fi  
    fi
}

    

  
  
	
 
