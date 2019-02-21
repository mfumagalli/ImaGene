
import os
import gzip

import numpy as np
import scipy.stats

import skimage.transform
from keras import models, layers, activations, optimizers, regularizers
from keras.utils import to_categorical, plot_model
from keras.models import load_model

import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import pymc3 # this will be removed
import pydot # optional

get_ipython().run_line_magic('run', '-i /rds/general/user/mfumagal/home/Software/ImaGene/ImaGene.py')

for i in 

myfile = ImaFile(simulations_folder='/home/mfumagal/Data/ImaGene.binary/Simulations1.Epoch3', nr_samples=128, model_name='Marth-3epoch-CEU')
