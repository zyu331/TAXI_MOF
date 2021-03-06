#!/bin/sh

#PBS -N SnapshotGCMC_AROFET
#PBS -A GT-dsholl3-joe
#PBS -o output
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -l walltime=2:00:00
#PBS -l mem=4gb
#PBS -S /bin/bash

cd $PBS_O_WORKDIR
#module purge
#module load intel-mpi/19.0.5
#module load python/3.7.4
module load anaconda3/2020.02
#module load openmpi/3.1.6
#module load gcc/10.1.0

time python <SnapshotGCMC.py> SnapshotGCMC.log
