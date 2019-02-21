#!/bin/bash

date

MODE=binary # or 'multi'

DIRSOFT=/home/mfumagal/Software
DIRDATA=/home/mfumagal/Data/ImaGene.$MODE

EACH=10
for PBS_ARRAY_INDEX in {0..29}
do
	A=$(($PBS_ARRAY_INDEX / $EACH))
	model=$(( A + 1 ))
	repetition=$(( $(( $PBS_ARRAY_INDEX - $(( A*EACH )) )) + 1 ))

	FNAME=$DIRDATA/Simulations$repetition.Epoch$model
	echo $FNAME
	mkdir -p $FNAME
	bash $DIRSOFT/ImaGene/simulate.sh $DIRSOFT/msms/lib/msms.jar $FNAME $model $MODE
done

date


