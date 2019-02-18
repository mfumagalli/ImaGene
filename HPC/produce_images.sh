#!/bin/bash

model=$1
repetition=$2

$IMAGENE=$HOME/Software/ImaGene

ipython $IMAGENE/HPC/produce_images.py $model $repetition

