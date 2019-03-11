
# simulate data for binary and multi-classification
bash generate_data.sh Binary

# training and testing for binary classification
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

# plot example images with different sorting
for s in 300 400
do
    for i in 0 1 2 3 4 5 6 7 8 9
    do
	python plot_examples.py $s $i
    done
done

# plot confusion matrices
python plot_binary.py

# plot training for one example
pyhton plot_train.py

# binary classification without dense layer
python train_binary_nodense.py > ~/Data/ImaGene/Logs/binary_nodense.txt

# binary classification without dense layer and 5x5 filter
python train_binary_5x5filter.py > ~/Data/ImaGene/Logs/binary_5x5filter.txt

# train and test for multiclassification for assessing model misrepresentation
for e in 3 2 1
do
    python train_multi.py $e > ~/Data/ImaGene/Logs/multi.$e.txt
done

# plot confusion matrices for multiclassification
python plot_multi.py

# simulate data for continuous parameter
bash generate_data.sh Continuous

# train and test for continuous parameter
python train_continuous.py 0 0 > ~/Data/ImaGene/Logs/continuous.0.0.txt
python train_continuous.py 0 0.5 > ~/Data/ImaGene/Logs/continuous.0.0.5.txt
python train_continuous.py 1 0 > ~/Data/ImaGene/Logs/continuous.1.0.txt

#python train_regression.py > ~/Data/ImaGene/Logs/regression.txt








