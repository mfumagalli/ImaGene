
bash generate_data.sh Binary

for e in 3 2 1
do
    for m in RowsCols Rows Cols None
    do
	    for s in 300 400 200
	    do
                python train_binary.py $e $m $s > ~/Data/ImaGene/Logs/binary.$e.$m.$s.txt
	done
    done
done

# generate one image of training loss!!! the other ones should use ALL the data, no validation set!

python analyse_binary.py

python train_binary_nodense.py > ~/Data/ImaGene/Logs/binary_nodense.txt

python train_binary_5x5filter.py > ~/Data/ImaGene/Logs/binary_5x5filter.txt

for e in 3 2 1
do
    python train_multi.py $e > ~/Data/ImaGene/Logs/multi.$e.txt
done




