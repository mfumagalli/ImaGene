
Pipeline on HPC system.

Perform simulations.
```
cp $HOME/Software/ImaGene/HPC/generate_data_binary_hpc.sh .
qsub generate_data_binary_hpc.sh
rm generate_data_binary_hpc.sh
```

```
cp $HOME/Software/ImaGene/HPC/generate_data_multi_hpc.sh .
qsub generate_data_multi_hpc.sh
rm generate_data_multi_hpc.sh
```


Produce images.
```
cp $HOME/Software/ImaGene/HPC/produce_images_hpc.sh .
qsub produce_images_hpc.sh
```




