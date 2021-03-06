name='CO2'
p=(10000 50000 100000 250000 500000 1000000)

j=0
for k in "${p[@]}"
	do
		mkdir $j
		cd $j
		sed "s,myP,$k,g" ../Snapshot.py > Snapshot.py
		sed "s,myP,$k,g" ../submit.sh > submit.sh
		cp ../data.* ../force_field* ../hMOFcifGen.py ../mof_taxi.py ../in.* .
		j=$((j+1))      

        qsub submit.sh
        cd ..
        
        done
      

