#!/bin/bash

date

DIRSOFT=/home/mfumagal/Software
DIRDATA=/home/mfumagal/Data/ImaGene

EACH=10
for PBS_ARRAY_INDEX in {0..29}
do
	A=$(($PBS_ARRAY_INDEX / $EACH))
	model=$(( A + 1 ))
	repetition=$(( $(( $PBS_ARRAY_INDEX - $(( A*EACH )) )) + 1 ))

	FNAME=$DIRDATA/Simulations$repetition.Epoch$model
	echo $FNAME
	mkdir -p $FNAME
	bash $DIRSOFT/ImaGene/Scripts/simulate.sh $DIRSOFT/msms/lib/msms.jar $FNAME $model
done

date


