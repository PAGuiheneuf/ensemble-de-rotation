#!/bin/bash
#
#PBS -l walltime=04:00:00
#PBS -N main
#PBS -j oe 
#PBS -l select=1:ncpus=32

module load mpt

cp $PBS_O_WORKDIR/rotationEnsemble .

export OMP_NUM_THREADS=32

omplace -nt 32 ./rotationEnsemble > log.txt

cp log.txt $PBS_O_WORKDIR
cp rotationEnsemble.txt $PBS_O_WORKDIR
