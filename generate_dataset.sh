
date

source $1

for (( INDEX=1; INDEX<=$NBATCH; INDEX++ ))
do
        FNAME=$DIRDATA/Simulations$INDEX
        echo $FNAME
        mkdir -p $FNAME

	for SEL in $SELRANGE
	do
		for TIME in $TIMERANGE
		do
    			java -jar $DIRMSMS -N $NREF -ms $NCHROMS $NREPL -t $THETA -r $RHO $LEN -Sp $SELPOS -SI $TIME 1 $FREQ -SAA $(($SEL*2)) -SAa $SEL -Saa 0 -Smark $DEMO -threads $NTHREADS | gzip > $FNAME/msms..$SEL..$TIME..txt.gz
		done
	done
done

date



