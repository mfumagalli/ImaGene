
# process simulations and produce images

# requirements for ImaGene

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

# for saving

import _pickle as pickle

# for input data

import sys

# 

model = str(sys.argv[1])
repetition = str(sys.argv[2])

# load ImaGene

get_ipython().run_line_magic('run', '-i /rds/general/user/mfumagal/home/Software/ImaGene/ImaGene.py')

# load file

dir_sim = '/rds/general/user/mfumagal/ephemeral/Data/ImaGene/Simulations' + repetition + '.Epoch' + model

myfile = ImaFile(simulations_folder=dir_sim, nr_samples=128, model_name=model)

mypop = myfile.read_simulations(parameter_name='selection_coeff_hetero')

print(mypop.summary())

# process

mypop.majorminor()

#mypop.filter_freq(0.01) # not done anymore

rnd_idx = get_index_random(mypop)

dir_name = '/rds/general/user/mfumagal/ephemeral/Data/ImaGene/Images' + repetition + '.Epoch' + model
if os.path.exists(dir_name) is False:
    os.mkdir(dir_name)

with open(dir_name + '/mypop','wb') as fp:
    pickle.dump(mypop, fp, protocol=4)
mypop.resize((128, 128))
mypop.convert()
mypop.subset(rnd_idx)
with open(dir_name + '/mypop_sortednone','wb') as fp:
    pickle.dump(mypop, fp, protocol=4)

with open(dir_name + '/mypop','rb') as fp:
    mypop = pickle.load(fp)
mypop.sort('rows_freq')
mypop.resize((128, 128))
mypop.convert()
mypop.subset(rnd_idx)
with open(dir_name + '/mypop_sortedrowsfreq','wb') as fp:
    pickle.dump(mypop, fp, protocol=4)

with open(dir_name + '/mypop','rb') as fp:
    mypop = pickle.load(fp)
mypop.sort('cols_freq')
mypop.resize((128, 128))
mypop.convert()
mypop.subset(rnd_idx)
with open(dir_name + '/mypop_sortedcolsfreq','wb') as fp:
    pickle.dump(mypop, fp, protocol=4)

with open(dir_name + '/mypop','rb') as fp:
    mypop = pickle.load(fp)
mypop.sort('rows_freq')
mypop.sort('cols_freq')
mypop.resize((128, 128))
mypop.convert()
mypop.subset(rnd_idx)
with open(dir_name + '/mypop_sortedrowscolsfreq','wb') as fp:
    pickle.dump(mypop, fp, protocol=4)

with open(dir_name + '/mypop','rb') as fp:
    mypop = pickle.load(fp)
mypop.sort('rows_dist')
mypop.sort('cols_dist')
mypop.resize((128, 128))
mypop.convert()
mypop.subset(rnd_idx)
with open(dir_name + '/mypop_sortedrowscolsdist','wb') as fp:
    pickle.dump(mypop, fp, protocol=4)


