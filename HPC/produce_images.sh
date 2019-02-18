#!/bin/bash

model=$1
repetition=$2

ipython $RDS/ephemeral/produce_images.py $model $repetition

