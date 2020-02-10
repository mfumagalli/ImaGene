#!/bin/bash

#PBS -l walltime=36:00:00
#PBS -l select=1:ncpus=8:mem=4gb
#PBS -J 0-29

date

MODE=binary # or 'multi'

DIRSOFT=/home/mfumagal/Software
DIRDATA=$EPHEMERAL/Data/ImaGene/$MODE

EACH=10
A=$(($PBS_ARRAY_INDEX / $EACH))
model=$(( A + 1 ))
repetition=$(( $(( $PBS_ARRAY_INDEX - $(( A*EACH )) )) + 1 ))

FNAME=$DIRDATA/Simulations$repetition.Epoch$model
echo $FNAME
mkdir -p $FNAME
bash $DIRSOFT/ImaGene/simulate.sh $DIRSOFT/msms/lib/msms.jar $FNAME $model $MODE

date


