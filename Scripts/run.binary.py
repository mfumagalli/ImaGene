
import os
import gzip

import numpy as np
import scipy.stats

import skimage.transform
from keras import models, layers, activations, optimizers, regularizers
from keras.utils import plot_model
from keras.models import load_model

import matplotlib.pyplot as plt
#from sklearn.metrics import confusion_matrix
#import pymc3 # this will be removed
import pydot # optional

import pathlib
import _pickle as pickle

exec(open('ImaGene.py').read())
#get_ipython().run_line_magic('run', '-i /home/mfumagal/Software/ImaGene/ImaGene.py')

for s in [100, 200, 300, 400]:

    for m in ['None', 'Rows', 'Cols', 'RowsCols']:

        for e in [1, 2, 3]:

            folder = '/home/mfumagal/Data/ImaGene/Binary/Results/Epoch' + str(e) + '/S' + str(s) + '/' + str(m)
            print(folder)
            pathlib.Path(folder).mkdir(parents=True, exist_ok=True) 
                       
            i = 0
            while i < 10:

                i += 1
                print(str(s) + str(m) + str(e) + str(i))

                myfile = ImaFile(simulations_folder='/home/mfumagal/Data/ImaGene/Binary/Simulations' + str(i) 
                         + '.Epoch' + str(e), nr_samples=128, model_name='Marth-' + str(e) + 'epoch-CEU')
                mypop = myfile.read_simulations(parameter_name='selection_coeff_hetero', max_nrepl=20)
    
                mypop.majorminor()
                mypop.filter_freq(0.01)
            
                if m == 'Rows':
                    mypop.sort('rows_freq')
                if m == 'Cols':
                    mypop.sort('cols_freq')
                if m == 'RowsCols':
                    mypop.sort('rows_freq')
                    mypop.sort('cols_freq')
                
                mypop.resize((128, 128))
                mypop.convert()
    
                mypop.classes = np.array([0,int(s)])
                classes_idx = get_index_classes(mypop.targets, mypop.classes)
                mypop.subset(classes_idx)
    
                rnd_idx = get_index_random(mypop)
                mypop.subset(rnd_idx)
    
                mypop.targets = to_binary(mypop.targets)
    
                if i == 1:
                    mynet = ImaNet(name='CPx2')
                    mynet.model = models.Sequential([
                        layers.Conv2D(filters=32, kernel_size=(3,3), strides=(1,1), activation='relu', kernel_regularizer=regularizers.l1_l2(l1=0.01, l2=0.01), padding='valid', input_shape=mypop.data.shape[1:4]),
                        layers.MaxPooling2D(pool_size=(2,2)),
                        #layers.Dropout(rate=0.5),
                        layers.Conv2D(filters=32, kernel_size=(3,3), strides=(1,1), activation='relu', kernel_regularizer=regularizers.l1_l2(l1=0.01, l2=0.01), padding='valid'),
                        layers.MaxPooling2D(pool_size=(2,2)),
                        #layers.Dropout(rate=0.5),
                        #layers.Conv2D(filters=32, kernel_size=(3,3), strides=(1,1), activation='relu', padding='valid'),
                        #layers.MaxPooling2D(pool_size=(2,2)),
                        #layers.Dropout(rate=0.5),
                        layers.Flatten(),
                        layers.Dense(units=64, activation='relu'),
                        layers.Dense(units=1, activation='sigmoid')])
                    mynet.model.compile(optimizer='rmsprop',
                        loss='binary_crossentropy',
                        metrics=['accuracy'])
                    mynet.plot_net(summary=True, file=folder + '/net.png')
                else:
                    model = load_model(folder + '/net.h5')
    
                if i < 10:
                    score = mynet.model.fit(mypop.data, mypop.targets, batch_size=32, epochs=1, verbose=0, validation_split=0.10)
                    mynet.update_scores(score)
                    mynet.model.save(folder + '/net.h5')
                else:
                    mynet.test = mynet.model.evaluate(mypop.data, mypop.targets, batch_size=None, verbose=0)
                    print(mynet.test)

            # save the latest data (testing data)
            with open(folder + '/mypop','wb') as fp:
                pickle.dump(mypop, fp)
            # save the latest network
            with open(folder + '/mynet','wb') as fp:
                pickle.dump(mynet, fp)
        
            del mypop
            del mynet






