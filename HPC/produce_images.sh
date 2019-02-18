#!/bin/bash

model=$1
repetition=$2

$DIRSOFT=$RDS/home/Software

ipython $DIRSOFT/ImaGene/HPC/produce_images.py $model $repetition

