#!/bin/sh

#PBS -N JN
#PBS -t 0-9
#PBS -A GT-dsholl3-joe
#PBS -o output
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20:00:00
#PBS -l mem=4gb
#PBS -S /bin/bash

cd $PBS_O_WORKDIR
tag=`printf "%s" $PBS_ARRAYID`
cd $tag

$RASPA_DIR/bin/simulate 
