#!/bin/sh

#PBS -N Snapshot_XXXXXX
#PBS -A GT-dsholl3-joe
#PBS -o output
#PBS -j oe
#PBS -l nodes=2:ppn=8
#PBS -l walltime=20:00:00
#PBS -l mem=4gb
#PBS -S /bin/bash

cd $PBS_O_WORKDIR
#module purge
#module load intel-mpi/19.0.5
#module load python/3.7.4
module load anaconda3/2020.02
#module load openmpi/3.1.6
#module load gcc/10.1.0

rm -r fig
mkdir fig
time python <Snapshot.py> Snapshot.log
