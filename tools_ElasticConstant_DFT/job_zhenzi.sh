#################
#PBS -N Elastic_XXX
#PBS -A GT-dsholl3-joe
#PBS -o output 
#PBS -j oe
#PBS -l nodes=1:ppn=24
#PBS -l walltime=300:00:00
#PBS -l mem=128gb
#PBS -S /bin/bash

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

module load intel/19.0.5 openmpi/3.1.6

mpirun -n 24 /storage/home/hcoda1/7/zyu331/p-dsholl3-0/bin/vasp.5.4.4.gam > out

###############
