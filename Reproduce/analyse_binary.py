
import numpy as np
import _pickle as pickle

import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import itertools

exec(open('/home/mfumagal/Software/ImaGene/ImaGene.py').read())

fig, ax = plt.subplots(3, 4, sharex='col', sharey='row')

for e in [3]:
    r = -1
    for s in [200,300,400]:
        r += 1
        c = -1
        classes = [0,s]
        classes = ['N','S']
        for m in ['None', 'Rows', 'Cols', 'RowsCols']:
            c += 1
            folder = '/home/mfumagal/Data/ImaGene/Binary/Results/Epoch' + str(e) + '/S' + str(s) + '/' + str(m)
            print(folder)
            mynet = load_imanet(folder + '/mynet')
            print(mynet.test)
            cm = confusion_matrix(mynet.values[0,:], mynet.values[1,:])
            accuracy = np.trace(cm) / float(np.sum(cm))
            cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
            ax[r,c].imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
            ax[r,c].set_xticks(np.arange(len(classes)))
            ax[r,c].set_yticks(np.arange(len(classes)))
            ax[r,c].set_xticklabels(classes)
            ax[r,c].set_yticklabels(classes)
            if r==0:
                if m == 'None':
                    ax[r,c].set_title('unsorted')
                if m == 'Rows':
                    ax[r,c].set_title('sorted by row')
                if m == 'Cols':
                    ax[r,c].set_title('sorted by column')
                if m == 'RowsCols':
                    ax[r,c].set_title('sorted by row and column')
            if r==2:
                ax[r,c].set_xlabel('Predicted')
            if c==0:
                ax[r,c].set_ylabel('True\nS = ' + str(s), rotation=0)
            thresh = cm.max() / 1.5
            for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
                ax[r,c].text(j, i, "{:0.3f}".format(cm[i, j]), horizontalalignment="center", color="white" if cm[i, j] > thresh else "black")

fig.show()


