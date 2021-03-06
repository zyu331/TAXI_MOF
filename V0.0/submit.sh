#!/bin/bash
#SBATCH --partition=hpg2-compute
#SBATCH --job-name=CO2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=2G
#SBATCH --time=96:00:00
#SBATCH --account=colina
#SBATCH --qos=colina

cd $SLURM_SUBMIT_DIR
module purge
module load intel/2018 openmpi/3.1.2
module load python/2.7.14
export LAMMPS_EXEC=/home/dylananstine/lammps/src/lmp_machine
time python <MCMD.py>MCMD.log
