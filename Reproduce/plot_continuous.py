
### --------------------------------- ###

# plot distributions

import numpy as np
import _pickle as pickle

import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix

import math

exec(open('/home/mfumagal/Software/ImaGene/ImaGene.py').read())

e = 3
folder = '/home/mfumagal/Data/ImaGene/Continuous/Results/Epoch' + str(e)

classes = [0,40,80,120,160,200,240,280,320,360,400]

fig, ax = plt.subplots(1, 3, sharex='col', sharey='row')

c = -1
for w in [0, 1]:
    for sd in [0, 0.5]:
        c += 1
        if (c<3):
            mynet = load_imanet(folder + '/Wiggle' + str(w) + '-Sd' + str(sd) + '/mynet')
            print(mynet.test)
            print(math.sqrt(np.average((mynet.values[0]-mynet.values[1])**2)))
            print(math.sqrt(np.average((mynet.values[0]-mynet.values[2])**2)))
            print(np.average((mynet.values[1]-mynet.values[0])/(mynet.values[0]+1)))
            print(np.average((mynet.values[2]-mynet.values[0])/(mynet.values[0]+1)))
            cm = confusion_matrix(mynet.values[0,:], mynet.values[1,:])
            accuracy = np.trace(cm) / float(np.sum(cm))
            cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
            ax[c].imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
            ax[c].set_xticks(np.arange(len(classes)))
            ax[c].set_yticks(np.arange(len(classes)))
            ax[c].set_xticklabels(classes)
            ax[c].set_yticklabels(classes)
            if c == 0:
                ax[c].set_title('Categorical')
            if c == 1:
                ax[c].set_title('Gaussian')
            if c == 2:
                ax[c].set_title('Categorical\nwith perturbation')
            ax[c].set_xlabel('Predicted')
            if c == 0:
                ax[c].set_ylabel('True')
            else:
                ax[c].set_ylabel('')

fig.show()


w = 0
sd = 0.5

mynet = load_imanet(folder + '/Wiggle' + str(w) + '-Sd' + str(sd) + '/mynet')
print(mynet.test)
print(math.sqrt(np.average((mynet.values[0]-mynet.values[1])**2)))
print(math.sqrt(np.average((mynet.values[0]-mynet.values[2])**2)))
print(np.average((mynet.values[2]-mynet.values[0])/(mynet.values[0])))

w = 1
sd = 0

mynet = load_imanet(folder + '/Wiggle' + str(w) + '-Sd' + str(sd) + '/mynet')
print(mynet.test)
print(math.sqrt(np.average((mynet.values[0]-mynet.values[1])**2)))
print(math.sqrt(np.average((mynet.values[0]-mynet.values[2])**2)))
print(np.average((mynet.values[2]-mynet.values[0])/(mynet.values[0])))

[1.7175514302883956, 0.3182543640897756]
64.99384360869304
55.23591758807737
[1.7587020619076088, 0.3129426433915212]
63.64113905390227
54.22718999510877
[1.886490782906587, 0.2416708229426434]
70.60691485811606
62.568657601569974

# best

from keras.models import load_model

import pymc3


w = 0
sd = 0.5

mygene = load_imagene(folder + '/Wiggle' + str(w) + '-Sd' + str(sd) + '/mygene')

model = load_model(folder + '/Wiggle' + str(w) + '-Sd' + str(sd) + '/model.h5')

mynet = load_imanet(folder + '/Wiggle' + str(w) + '-Sd' + str(sd) + '/mynet')

plt.figure()

plt.subplot(121)

ii = 3

ind = np.where( ((mynet.values[0]-mynet.values[2])**2<1) & (mynet.values[0]>100) & (mynet.values[0]<130))
probs = model.predict(mygene.data[ind[0][ii]].reshape(1,128,128,1), batch_size=None)[0]
samples_distr = np.random.choice(mygene.classes, size = 100000, replace = True, p = probs)
HPD = pymc3.stats.hpd(samples_distr, alpha = 0.05)
BF = (1 - probs[0]) / probs[0]
MAP = mygene.classes[np.argmax(probs)]
MLE = np.average(mygene.classes, weights = probs)
print(BF)
tick_marks = mygene.classes
cen_tick = mygene.classes
plt.hist(samples_distr, color='#a6bddb', bins=len(mygene.classes), density=True)
plt.xlim([mygene.classes.min(), mygene.classes.max()])
plt.xticks(cen_tick, cen_tick, rotation=45, fontsize=10)
plt.yticks(fontsize=10)
plt.ylabel('Density', fontsize=12)
plt.xlabel('Parameter', fontsize=12)
plt.title('Sampled posterior distribution')
plt.grid(True)
plt.axvline(MLE, label='mean ('+str(round(MLE,2))+')', color='r', linestyle='--')
plt.axvline(MAP, label='MAP ('+str(MAP)+')', color='b', linestyle='--')
plt.axhline(y=0.0001, xmin=HPD[0]/np.max(mygene.classes), xmax=HPD[1]/np.max(mygene.classes), c='black', label='95% HPD\nInterval: [{}, {}]'.format(HPD[0],HPD[1]))
plt.legend()

plt.subplot(122)

ii = 0

ind = np.where( ((mynet.values[0]-mynet.values[2])**2<1) & (mynet.values[0]>300) & (mynet.values[0]<330))
probs = model.predict(mygene.data[ind[0][ii]].reshape(1,128,128,1), batch_size=None)[0]
samples_distr = np.random.choice(mygene.classes, size = 100000, replace = True, p = probs)
HPD = pymc3.stats.hpd(samples_distr, alpha = 0.05)
BF = (1 - probs[0]) / probs[0]
MAP = mygene.classes[np.argmax(probs)]
MLE = np.average(mygene.classes, weights = probs)
print(BF)
tick_marks = mygene.classes
cen_tick = mygene.classes
plt.hist(samples_distr, color='#a6bddb', bins=len(mygene.classes), density=True)
plt.xlim([mygene.classes.min(), mygene.classes.max()])
plt.xticks(cen_tick, cen_tick, rotation=45, fontsize=10)
plt.yticks(fontsize=10)
plt.ylabel('Density', fontsize=12)
plt.xlabel('Parameter', fontsize=12)
plt.title('Sampled posterior distribution')
plt.grid(True)
plt.axvline(MLE, label='mean ('+str(round(MLE,2))+')', color='r', linestyle='--')
plt.axvline(MAP, label='MAP ('+str(MAP)+')', color='b', linestyle='--')
plt.axhline(y=0.0001, xmin=HPD[0]/np.max(mygene.classes), xmax=HPD[1]/np.max(mygene.classes), c='black', label='95% HPD\nInterval: [{}, {}]'.format(HPD[0],HPD[1]))
plt.legend()

plt.show()



19.86593997388371
87235.51817515353









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
    ax[c].set_title('Tested on ' + str(e) + '-epoch model')
    ax[c].set_xlabel('Predicted')
    if c == 0:
        ax[c].set_ylabel('True\nTrained on 3-epoch model')
    else:
        ax[c].set_ylabel('')
    thresh = cm.max() / 1.5
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        ax[c].text(j, i, "{:0.3f}".format(cm[i, j]), horizontalalignment="center", color="white" if cm[i, j] > thresh else "black")

fig.show()

plt.savefig(fname = '/home/mfumagal/Data/ImaGene/Multi/Results/Figure_demo.pdf')

### --------------------------------- ###






