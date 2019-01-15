
### Configuration for simulations 

# bash script.sh msms.lib folder model

### ---------------------------------------------------------------

mkdir -p $2

## 1) DEMOGRAPHIC MODEL

# You need to specify the demographic model to be used in your simulations. The model should be scaled used a reference effective population size.
# As an illustration, the following parameters for a demographic model are taken from Marth et al. Genetics 2004 (table 1, European data, 1,2, or 3 epoch model):
# three epoch: N3=10,000,N2=2,000,N1=20,000,T2=5,00,T1=3,000
# one epoch: N1=10,000
# two epoch: N2=10,000, N1=140,000, T1=2000

NREF=10000 # reference effective population size

MARTH1='' # Marth 1-epoch for CEU
MARTH2='-eN 0.05 1 -eN 0 14' # Marth 2-epoch for CEU
MARTH3='-eN 0.0875 1 -eN 0.075 0.2 -eN 0 2' # Marth 3-epoch for CEU

# Once you decide (or implement) the demographic model to be used, change the following variable accordingly.
# As an illustration, here we used a 3-epoch model for Europeans

if [ "$3" -eq "1" ]; then
	DEMO=$MARTH1;
fi
if [ "$3" -eq "2" ]; then
	DEMO=$MARTH2;
fi
if [ "$3" -eq "3" ]; then
	DEMO=$MARTH3;
fi


## 2) LOCUS AND SAMPLE SIZE

# Then you need to define parameters on the genomic locus to estimate. These parameters are the length (in bp), mutation and recombination rate. You should need to specify how many samples (e.g. haplotypes) you wish to extract from the simulations. Be aware that mutation and recombination rates are scaled in 4*Ne*length, and therefore you need to use the sample NREF specified in the step above.

LEN=100000 # length of the locus in bp
THETA=60 # mutation rate in 4*Ne*LEN scale; 60 corresponds to 1.5e-8 for Ne=10,000 and 100,000 bp length
RHO=40 # recombination rate (rho); 40 corresponds to 1e-8 for Ne=10,000 and 100,000 bp length

NCHROMS=128 # number of haplotypes (chromosomes) to extract: 198 matches the number of  unrelated CEU samples in 1000 Genomes Project data

## 3) SELECTION

# Finally, you need to specify parameters for the positive selection event. These include the time of onset, the position of selected allele, initial frequency, and the range for the selection coefficient to be estimated. Be aware that selection time is scaled in 4*Ne generations while the selection coefficient is in 2*Ne units. The selection coefficient is at the allele level and we assume an additive effect.

SELPOS=`bc <<< 'scale=2; 1/2'` # relative position of selected allele; the example here indicates that the selected allele sits in the middle of the locus

FREQ=`bc <<< 'scale=6; 1/20000'` # frequency of selected allele at start of selection; here de novo, 1/2N

SELRANGE=`seq 0 1 799` # range and step for the selection coefficient to be estimated in 2*Ne units;

#NREPL=125 # this is the number of replicates (simulations) per value of selection coefficient to be estimated
NREPL=10

SELTIME=`bc <<< 'scale=4; 600/40000'` # 15kya

# time for the start of selection in 4*Nref generations; e.g. 800/40000 is at 20kya, with Ne=10k and 25 years as gen time.
for SEL in $SELRANGE
do
    java -jar $1 -N $NREF -ms $NCHROMS $NREPL -t $THETA -r $RHO $LEN -Sp $SELPOS -SI $SELTIME 1 $FREQ -SAA $(($SEL*2)) -SAa $SEL -Saa 0 -Smark $DEMO -thread 4 | gzip > $2/msms..$SEL..$SELTIME..txt.gz
done



