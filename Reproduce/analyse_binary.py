
import os
import gzip
import _pickle as pickle

import numpy as np
import scipy.stats

import skimage.transform
from keras import models, layers, activations, optimizers, regularizers
from keras.utils import plot_model
from keras.models import load_model
from keras import backend as K

import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import pymc3 # this will be removed
import pydot # optional

exec(open('/home/mfumagal/Software/ImaGene/ImaGene.py').read())

import sys

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

e = str(sys.argv[1]) # epoch
m = str(sys.argv[2]) # sorting
s = str(sys.argv[3]) # sel coeff

folder = '/home/mfumagal/Data/ImaGene/Binary/Results/Epoch' + str(e) + '/S' + str(s) + '/' + str(m)

model = load_model(folder + '/model.h5')

mygene = load_imagene(folder + '/mygene')

mynet = load_imanet(folder + '/mynet')







