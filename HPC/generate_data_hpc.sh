#!/bin/bash

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=4:mem=1gb
#PBS -J 0-29

date

DIRSOFT=$RDS/home/Software
DIRDATA=$RDS/ephemeral/Data/ImaGene

EACH=10
A=$(($PBS_ARRAY_INDEX / $EACH))
model=$(( A + 1 ))
repetition=$(( $(( $PBS_ARRAY_INDEX - $(( A*EACH )) )) + 1 ))

FNAME=$DIRDATA/Simulations$repetition.Epoch$model
echo $FNAME
mkdir -p $FNAME
bash $DIRSOFT/ImaGene/Scripts/simulate.sh $DIRSOFT/msms/lib/msms.jar $FNAME $model

date


