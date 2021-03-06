#!/bin/sh

#PBS -N CAU-10-1bar
#PBS -A GT-dsholl3-joe
#PBS -o output
#PBS -j oe
#PBS -l nodes=2:ppn=16
#PBS -l walltime=100:00:00
#PBS -l mem=4gb
#PBS -S /bin/bash

cd $PBS_O_WORKDIR
#module purge
#module load intel-mpi/19.0.5
#module load python/3.7.4
module load anaconda2/2019.10
#module load openmpi/3.1.6
#module load gcc/10.1.0

time python <MCMD.py> MCMD.log
