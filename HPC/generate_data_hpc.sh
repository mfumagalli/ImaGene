#!/bin/bash

#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=4:mem=1gb
#PBS -J 0-29

date

IMAGENE=$HOME/Software/ImaGene
MSMS=$HOME/Software/msms

DATA=$EPHEMERAL/Data/ImaGene

EACH=10
A=$(($PBS_ARRAY_INDEX / $EACH))
model=$(( A + 1 ))
repetition=$(( $(( $PBS_ARRAY_INDEX - $(( A*EACH )) )) + 1 ))

FNAME=$DATA/Simulations$repetition.Epoch$model
echo $FNAME
mkdir -p $FNAME
bash $IMAGENE/simulate.sh $MSMS/lib/msms.jar $FNAME $model

date


