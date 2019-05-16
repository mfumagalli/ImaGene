
# simulate binary/multiclass (option 'Binary') or (semi)-continuous (option 'Continuous')
MODE=$1 # Binary or Continuous

# path to software (ImaGene and msms) and where to store simulations
DIRSOFT=/home/mfumagal/Software
DIRDATA=/home/mfumagal/Data/ImaGene/$MODE

date

# number of batches of simulations
EACH=10
# each batch is simulated 3 times, for each demographic model
for INDEX in {0..29}
do
	A=$(($INDEX / $EACH))
	model=$(( A + 1 ))
	repetition=$(( $(( $INDEX - $(( A*EACH )) )) + 1 ))

	FNAME=$DIRDATA/Simulations$repetition.Epoch$model
	echo $FNAME
	mkdir -p $FNAME
	bash $DIRSOFT/ImaGene/Reproduce/simulate.sh $DIRSOFT/msms/lib/msms.jar $FNAME $model $MODE
done

date


