#This file contains the defaults for compiling FLEUR on known machines

#please add further definitions here and also add code below!
read -r -d '' known_machines << EOM
   IFF         -- cluster @ PGI of FZJ
   JURECA      -- @JSC
EOM


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

    

  
  
	
 
