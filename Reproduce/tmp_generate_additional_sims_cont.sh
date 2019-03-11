#!/bin/bash

date

MODE=Continuous

DIRSOFT=/home/mfumagal/Software
DIRDATA=/home/mfumagal/Data/ImaGene/$MODE

model=3

for repetition in {11..50}
do

	FNAME=$DIRDATA/Simulations$repetition.Epoch$model
	echo $FNAME
	mkdir -p $FNAME
	bash $DIRSOFT/ImaGene/Reproduce/simulate.sh $DIRSOFT/msms/lib/msms.jar $FNAME $model $MODE
done

date


