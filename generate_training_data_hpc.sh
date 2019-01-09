
for repetition in 1
do
        for model in 1 2 3
        do
                FNAME=$RDS/ephemeral/Data/ImaGene/Simulations$repetition.Epoch$model
                mkdir -p $FNAME
                echo $FNAME
                bash $RDS/home/Software/ImaGene/Scripts/simulations_denovo_continuous.sh $RDS/home/Software/msms/lib/msms.jar $FNAME $model
        done
done


