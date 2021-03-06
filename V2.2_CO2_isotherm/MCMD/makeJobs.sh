name='CO2'
p=(10000 100000 250000 500000 1000000 1500000)

j=0
for k in "${p[@]}"
	do
		mkdir $j
		cd $j
		sed "s,myP,$k,g" ../MCMD.py > MCMD.py
		sed "s,myP,$k,g" ../submit.sh > submit.sh
		cp ../data.* ../force_field* ../pseudo_atoms.def ../mof_taxi.py ../in.* ../CO2.def ../clean.sh ../*.cif .
		j=$((j+1))      

        qsub submit.sh
        cd ..
        
        done
      

