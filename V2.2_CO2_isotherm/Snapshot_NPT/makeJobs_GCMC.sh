name='CO2'
p=(10000 50000 100000 250000 500000 1000000)

j=0
for k in "${p[@]}"
do
	cd $j
	mkdir MD
	mv * MD
	cp -r MD/cifs .
	mkdir GCMC
	cp ../force* ../pse* GCMC
	sed "s,myP,$k,g" ../simulation.input > GCMC/simulation.input 
	sed "s,JN,$k,g" ../JobArray_Ad.sh > GCMC/JobArray_Ad.sh
	cd GCMC
	i=0
	while [ $i -lt 9 ]; 
		do

		mkdir $i
		sed "s,mya,$i,g" simulation.input > $i/simulation.input
		cp ../cifs/${i}.cif force* pse* $i
		i=$((i+1))

		done
	qsub JobArray_Ad.sh
	cd ..

	j=$((j+1))   
        cd ..
done
      

