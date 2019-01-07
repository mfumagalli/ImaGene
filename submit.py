
# generate data

import os
import gzip

import numpy as np
import scipy.stats

import skimage.transform
from keras import models, layers, optimizers, regularizers
from keras.utils import to_categorical, plot_model

import matplotlib.pyplot as plt
import pymc3
import pydot

# 

import _pickle as pickle

#
get_ipython().run_line_magic('run', '-i /home/mfumagal/Software/ImaGene/ImaGene.py')

myfile = ImaFile(simulations_folder='/home/mfumagal/Data/ImaGene/Simulations1.Epoch1', nr_samples=128, model_name='Marth-1epoch-CEU')

mypop = myfile.read_simulations(parameter_name='selection_coeff_hetero', max_nrepl=10)

print(mypop.summary())

mypop.majorminor()

mypop.filter_freq(0.01)

shuffle_index = np.random.permutation(len(mypop.data))

os.mkdir('/home/mfumagal/Data/ImaGene/Sorting_effect_Epoch1')

print('saving')
with open('/home/mfumagal/Data/ImaGene/Sorting_effect_Epoch1/mypop','wb') as fp:
    pickle.dump(mypop, fp)
print('saved')

mypop.resize((128, 128))
mypop.convert()
mypop.shuffle(shuffle_index)
print('saving')
with open('/home/mfumagal/Data/ImaGene/Sorting_effect_Epoch1/mypop_sortednone','wb') as fp:
    pickle.dump(mypop, fp)
print('saved')

with open('/home/mfumagal/Data/ImaGene/Sorting_effect_Epoch1/mypop','rb') as fp:
    mypop = pickle.load(fp)
mypop.sort('rows_freq')
mypop.resize((128, 128))
mypop.convert()
mypop.shuffle(shuffle_index)
print('saving')
with open('/home/mfumagal/Data/ImaGene/Sorting_effect_Epoch1/mypop_sortedrowsfreq','wb') as fp:
    pickle.dump(mypop, fp)
print('saved')

with open('/home/mfumagal/Data/ImaGene/Sorting_effect_Epoch1/mypop','rb') as fp:
    mypop = pickle.load(fp)
mypop.sort('cols_freq')
mypop.resize((128, 128))
mypop.convert()
mypop.shuffle(shuffle_index)
print('saving')
with open('/home/mfumagal/Data/ImaGene/Sorting_effect_Epoch1/mypop_sortedcolsfreq','wb') as fp:
    pickle.dump(mypop, fp)
print('saved')

with open('/home/mfumagal/Data/ImaGene/Sorting_effect_Epoch1/mypop','rb') as fp:
    mypop = pickle.load(fp)
mypop.sort('rows_freq')
mypop.sort('cols_freq')
mypop.resize((128, 128))
mypop.convert()
mypop.shuffle(shuffle_index)
print('saving')
with open('/home/mfumagal/Data/ImaGene/Sorting_effect_Epoch1/mypop_sortedrowscolsfreq','wb') as fp:
    pickle.dump(mypop, fp)
print('saved')

with open('/home/mfumagal/Data/ImaGene/Sorting_effect_Epoch1/mypop','rb') as fp:
    mypop = pickle.load(fp)
mypop.sort('rows_distance_top')
mypop.sort('cols_distance_top')
mypop.resize((128, 128))
mypop.convert()
mypop.shuffle(shuffle_index)
print('saving')
with open('/home/mfumagal/Data/ImaGene/Sorting_effect_Epoch1/mypop_sortedrowscolsdist','wb') as fp:
    pickle.dump(mypop, fp)
print('saved')


