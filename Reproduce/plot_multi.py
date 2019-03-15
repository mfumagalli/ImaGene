
### --------------------------------- ###

# plot confusion matrices for multiclass classification

import numpy as np
import _pickle as pickle

import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import itertools

exec(open('/home/mfumagal/Software/ImaGene/ImaGene.py').read())

fig, ax = plt.subplots(1, 2, sharex='col', sharey='row')

classes = [0,200,400]
m = 'RowsCols'

c = -1
for e in [3,1]:

    c += 1
    folder = '/home/mfumagal/Data/ImaGene/Multi/Results/Epoch' + str(e) 
    print(folder)
    mynet = load_imanet(folder + '/mynet')
    print(mynet.test)
    cm = confusion_matrix(mynet.values[0,:], mynet.values[1,:])
    accuracy = np.trace(cm) / float(np.sum(cm))
    cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    ax[c].imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
    ax[c].set_xticks(np.arange(len(classes)))
    ax[c].set_yticks(np.arange(len(classes)))
    ax[c].set_xticklabels(classes)
    ax[c].set_yticklabels(classes)
    ax[c].set_title('Trained on ' + str(e) + '-epoch model')
    ax[c].set_xlabel('Predicted')
    if c == 0:
        ax[c].set_ylabel('Tested on 3-epoch model\nTrue')
    else:
        ax[c].set_ylabel('')
    thresh = cm.max() / 1.5
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        ax[c].text(j, i, "{:0.3f}".format(cm[i, j]), horizontalalignment="center", color="white" if cm[i, j] > thresh else "black")

fig.show()

plt.savefig(fname = '/home/mfumagal/Data/ImaGene/Multi/Results/Figure_demo.pdf')

### --------------------------------- ###






