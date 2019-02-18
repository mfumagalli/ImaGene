!/bin/bash

#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=1:mem=64gb
#PBS -J 0-29

date

EACH=10
A=$(($PBS_ARRAY_INDEX / $EACH))
model=$(( A + 1 ))
repetition=$(( $(( $PBS_ARRAY_INDEX - $(( A*EACH )) )) + 1 ))

echo $model $repetition

module load anaconda3/personal

bash $RDS/ephemeral/produce_images.sh $model $repetition

date

