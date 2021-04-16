
p=(ethene ethane propene propane methane xe kr)


for k in "${p[@]}"
	do
		cd $k
       	 qsub submit.sh
        	cd ..
        
        done
      

