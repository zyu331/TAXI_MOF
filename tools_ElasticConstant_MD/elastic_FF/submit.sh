#!/bin/sh

#PBS -N Elasicity_UIO-66
#PBS -A GT-dsholl3-joe
#PBS -o output
#PBS -j oe
#PBS -l nodes=2:ppn=4
#PBS -l walltime=20:00:00
#PBS -l mem=4gb
#PBS -S /bin/bash

cd $PBS_O_WORKDIR

rm out* log*
module load anaconda3
python elasticPreparation.py

mpirun -np 8 $LAMMPS_EXEC < in.elastic > out
