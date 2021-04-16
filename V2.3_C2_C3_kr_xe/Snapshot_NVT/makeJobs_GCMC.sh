N='myAd'

mkdir MD
mv * MD
cp -r MD/cifs .
mkdir GCMC

cp ../force* ../pse* GCMC
cp  MD/*.def GCMC
sed "s,Adorbate,$N,g" ../simulation.input > GCMC/simulation.input
sed "s,JN,$N,g" ../JobArray_Ad.sh > GCMC/JobArray_Ad.sh
cd GCMC

i=0
mkdir $i
sed "s,mya,relaxed,g" simulation.input > $i/simulation.input
cp ../cifs/relaxed.cif force* pse* $i
cp *.def $i
i=$((i+1))
	
while [ $i -lt 10 ]; 
	do

	mkdir $i
	sed "s,mya,$i,g" simulation.input > $i/simulation.input
	cp ../cifs/${i}.cif force* pse* $i
	cp *.def $i
	i=$((i+1))

	done
qsub JobArray_Ad.sh


  


      

