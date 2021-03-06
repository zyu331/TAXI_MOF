name='CO2'
p=(10000 50000 100000 250000 500000 1000000)

i=0
while [ $i -lt 9 ]; do
	mkdir $i
	cd $i
	j=0
	for k in "${p[@]}"
	do
		mkdir $j
		sed "s,mya,$i,g" ../simulation.input > $j/simulation.input
		sed -i "s,myP,$k,g" $j/simulation.input
		cp ../cifs/${i}.cif ../force_* ../pse* $j/
		j=$((j+1))
        done
        
        sed "s,JN,$name,g" ../JobArray_Ad.sh > JobArray_Ad.sh
        qsub JobArray_Ad.sh
        cd ..
        i=$((i+1))
        
done

