
# simulate selection under 3 different models

for repetition in 1
do
	for model in 1 2 3 
	do
    		FNAME=/home/mfumagal/Data/ImaGene/Simulations$repetition.Epoch$model
    		mkdir -p $FNAME
    		echo $FNAME
    		bash Scripts/simulations.sh /home/mfumagal/Software/msms/lib/msms.jar $FNAME $model
	done
done


# de novo, 15kya

for repetition in 2
do
        for model in 1 2 3
        do
                FNAME=/home/mfumagal/Data/ImaGene/Simulations$repetition.Epoch$model
                mkdir -p $FNAME
                echo $FNAME
                bash Scripts/simulations_denovo.sh /home/mfumagal/Software/msms/lib/msms.jar $FNAME $model
        done
done

# de novo continuous

for repetition in 3
do
        for model in 1 2 3
        do
                FNAME=/home/mfumagal/Data/ImaGene/Simulations$repetition.Epoch$model
                mkdir -p $FNAME
                echo $FNAME
                bash Scripts/simulations_denovo_continuous.sh /home/mfumagal/Software/msms/lib/msms.jar $FNAME $model
        done
done


