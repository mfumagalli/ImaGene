
### --------------------------------- ###

# plot training for binary classification

import numpy as np
import _pickle as pickle

import matplotlib.pyplot as plt

exec(open('/home/mfumagal/Software/ImaGene/ImaGene.py').read())

e = 3
s = 300
m = 'RowsCols'

folder = '/home/mfumagal/Data/ImaGene/Binary/Results/Epoch' + str(e) + '/S' + str(s) + '/' + str(m)
mynet = load_imanet(folder + '/mynet')

mynet.plot_train()


### --------------------------------- ###






