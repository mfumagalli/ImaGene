#!/bin/bash

#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=2:mem=256gb

DIRSOFT=$HOME/Software/ImaGene

module load anaconda3/personal

date

python train_binary_hpc.py > log_binary

date

