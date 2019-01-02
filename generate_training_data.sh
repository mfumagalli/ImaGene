
# simulate selection at 30kya under 3 different models

for model in 1 2 3 
do
    FNAME=/home/mfumagal/Data/ImaGene/Simulations_Epoch$model
    mkdir -p $FNAME
    echo $FNAME
    bash Scripts/simulations.sh /home/mfumagal/Software/msms/lib/msms.jar $FNAME $model
done

