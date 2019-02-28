#!/bin/bash

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=64gb
#PBS -J 0-2

date

EPOCH=$(($PBS_ARRAY_INDEX + 1))

IMAGENE=$HOME/Software/ImaGene

module load anaconda3/personal

python $IMAGENE/HPC/train_binary_hpc.py $EPOCH > log_binary

date

