
# plot one example of sorted images

import os
import gzip
import skimage.transform

import numpy as np
import matplotlib.pyplot as plt

import sys

exec(open('/home/mfumagal/Software/ImaGene/ImaGene.py').read())

e = 3
i = 1
s = str(sys.argv[1])
index = str(sys.argv[2])

fig, ax = plt.subplots(2, 2, sharex='col', sharey='row')

for r in range(2):

    for c in range(2):

        myfile = ImaFile(simulations_folder='/home/mfumagal/Data/ImaGene/Binary/Simulations' + str(i) + '.Epoch' + str(e), nr_samples=128, model_name='Marth-' + str(e) + 'epoch-CEU')
        mygene = myfile.read_simulations(parameter_name='selection_coeff_hetero', max_nrepl=10)

        mygene.majorminor()
        mygene.filter_freq(0.01)

        if (r==0) & (c==0):
            m = 'None'
        if (r==0) & (c==1):
            m = 'Rows'
        if (r==1) & (c==0):
            m = 'Cols'
        if (r==1) & (c==1):
            m = 'RowsCols'

        if (m =='Rows') | (m == 'RowsCols'):
            mygene.sort('rows_freq')
        if (m =='Cols') | (m == 'RowsCols'):
            mygene.sort('cols_freq')

        mygene.resize((128, 128))
        mygene.convert()

        ax[r,c].imshow(mygene.data[int(index)][:,:,0], cmap='gray')
        if r == 1:
            ax[r,c].set_xlabel('Position')
        if c == 0:
            ax[r,c].set_ylabel('Haplotype', rotation=0)

        if m == 'None':
            ax[r,c].set_title('unsorted')
        if m == 'Rows':
            ax[r,c].set_title('sorted by row')
        if m == 'Cols':
            ax[r,c].set_title('sorted by column')
        if m == 'RowsCols':
            ax[r,c].set_title('sorted by row and column')


fig.show()

fname = '/home/mfumagal/Data/ImaGene/Binary/Results/Epoch' + str(e) + '/S' + str(s) + '/example' + str(index) + '.pdf'
print(fname)

plt.savefig(fname=fname) 



