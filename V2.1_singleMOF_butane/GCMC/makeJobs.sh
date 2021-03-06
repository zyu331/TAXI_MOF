j='butane'
mkdir Henry
mkdir Adsorption
i=0
while [ $i -lt 10 ]; do
	mkdir Henry/$i
	sed "s,mya,$i,g" simulation.input_Henry > Henry/$i/simulation.input
	cp ../cifs/${i}.cif force_* pse* Henry/$i

	mkdir Adsorption/$i
	sed "s,mya,$i,g" simulation.input > Adsorption/$i/simulation.input
	cp ../cifs/${i}.cif force_* pse* Adsorption/$i
        i=$((i+1))
done

sed "s,JN,$j,g" JobArray_Henry.sh > Henry/JobArray.sh

sed "s,JN,$j,g" JobArray_Ad.sh > Adsorption/JobArray.sh

