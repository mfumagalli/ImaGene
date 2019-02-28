
bash generate_data.sh Binary

for e in 3 2 1
do
    for m in RowsCols Rows Cols None
    do
	    for s in 300 400 200
	    do
                python train_binary.py $e $m $s > Logs/binary.$e.$m.$s.txt
	done
    done
done




