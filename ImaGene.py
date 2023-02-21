
### ------------- utilities --------------------


def extract_msms_parameter(line, option, position=0):
    """
    Extract simulation parameters from first line of msms file

    Keyword Arguments:
        line (string) -- first line of (gzipped) msms file in bytes format 
        option (string) -- switch of msms to match
        position (int) -- i-th value to be taken after the switch

    Return:
        parameter (string)
    """

    return line.partition(option)[2].split()[position]


def get_index_classes(targets, classes):
    """
    Get index array for targets corresponding to selected classes

    Keyword Arguments:
        targets (array) -- target feature in ImaGene object
        classes (array) -- classes to select from targets

    Return:
        index (array)
    """
    index = []
    for counter,value in enumerate(classes):
        index = np.concatenate([index, np.where(targets==value)[0]])
    return np.asarray(index, dtype='int')


def get_index_random(genes=[], length=0):
    """
    Get random index array

    Keyword Arguments:
        length (int) -- length random index array
        genes (object) -- ImaGene object

    Return:
        index (array)
    """
    if length == 0:
        if len(genes.data) == 0:
            print('Either length or genes must be provided.')
        else:
            length = len(genes.data)

    return np.random.permutation(length)


def calculate_allele_frequency(genes, position):
    """
    ...
    """
    return [np.where(genes.data[i][:,np.where(genes.positions[i]==position)[0][0],0]==255,1,0).sum() for i in range(len(genes.data))]


def to_binary(targets):
    return np.asarray(np.where(targets == targets.min(), 0, 1).astype('float32'))


def to_categorical(targets, wiggle=0, sd=0):
    classes = np.unique(targets)
    nr_classes = len(classes)
    results = np.zeros((len(targets), len(classes)), dtype='float32')
    for counter, value in enumerate(targets):
        index = np.where(classes == value)[0]
        # add wiggle (if any)
        if wiggle > 0:
            index += np.random.randint(low=-wiggle, high=wiggle+1)
            if index < 0:
                index = 0
            elif index >= results.shape[1]:
                index = results.shape[1] - 1
        results[counter, index] = 1.
        # add sd (if any)
        if sd > 0:
            probs = scipy.stats.norm.pdf(range(nr_classes), loc=index, scale=sd)
            results[counter, ] = probs / probs.sum()
            del probs
    return results

def load_imagene(file):
    """
    Load ImaGene object
    """
    with open(file, 'rb') as fp:
        gene = pickle.load(fp)
    return gene

def load_imanet(file):
    """
    Load ImaNet object
    """
    with open(file, 'rb') as fp:
        net = pickle.load(fp)
    return net

def plot_scores(model, gene, classes, H0_class=0):
    """
    Plot scores of a predicted image as posterior distribution
    """
    probs = model.predict(gene.data, batch_size=None)[0]
    # Monte Carlo sampling
    samples_distr = np.random.choice(classes, size = 100000, replace = True, p = probs)
    # summary statistics and metrics of confidence
    # HPD = pymc3.stats.hpd(samples_distr, credible_interval = 0.95)
    HPD = az.hdi(samples_distr, credible_interval = 0.95)
    BF = (1 - probs[H0_class]) / probs[H0_class]
    MAP = classes[np.argmax(probs)]
    MLE = np.average(classes, weights = probs)
    
    # plot
    tick_marks = classes
    cen_tick = classes
    plt.hist(samples_distr, color='#a6bddb', bins=len(classes), density=True)
    plt.xlim([classes.min(), classes.max()])
    plt.xticks(cen_tick, cen_tick, rotation=45, fontsize=10)
    plt.yticks(fontsize=10)
    plt.ylabel('Density', fontsize=12)
    plt.xlabel('Parameter', fontsize=12)
    plt.title('Sampled posterior distribution')
    plt.grid(True)
    plt.axvline(MLE, label='mean ('+str(round(MLE,2))+')', color='r', linestyle='--')
    plt.axvline(MAP, label='MAP ('+str(MAP)+')', color='b', linestyle='--')
    plt.axhline(y=0.0001, xmin=HPD[0]/np.max(classes), xmax=HPD[1]/np.max(classes), c='black', label='95% HPD\nInterval: [{}, {}]'.format(HPD[0],HPD[1]))
    plt.legend()

    return (MAP, MLE, HPD, BF)






### -------- objects ------------------

class ImaFile:
    """
    Parser for real data and simulations
    """
    def __init__(self, nr_samples, simulations_folder=None, VCF_file_name=None, model_name='N/A'):
        self.simulations_folder = simulations_folder
        self.nr_samples = nr_samples
        self.VCF_file_name = VCF_file_name
        self.model_name = model_name
        return None

    def extract_description(self, file_name, first_line):
        """
        Read first line of simulations, extract all metadata and store it in a dictionary

        Keyword Arguments:
            file_name (string) -- name of simulation file
            first_line (string) -- first line of gzipped msms file
            model_name (string) -- name of demographic model

        Return:
            description (string)
        """

        desc = {'name':file_name}

        # Extracting parameters
        desc.update({'Nref':int(extract_msms_parameter(first_line, '-N '))})
        desc.update({'nr_chroms':int(extract_msms_parameter(first_line, '-N ', 1))})
        desc.update({'nr_replicates':int(extract_msms_parameter(first_line, '-N ', 2))})

        desc.update({'mutation_rate':float(extract_msms_parameter(first_line, '-t '))})
        desc.update({'recombination_rate':float(extract_msms_parameter(first_line, '-r '))})
        desc.update({'recombination_rate_nr_sites':int(extract_msms_parameter(first_line, '-r ', 1))})

        desc.update({'selection_position':float(extract_msms_parameter(first_line, '-Sp '))})
        desc.update({'selection_start_time':float(extract_msms_parameter(first_line, '-SI '))})
        desc.update({'selection_start_frequency':float(extract_msms_parameter(first_line, '-SI ', 2))})
    
        desc.update({'selection_coeff_HOMO':int(extract_msms_parameter(first_line, '-SAA '))})
        desc.update({'selection_coeff_hetero':int(extract_msms_parameter(first_line, '-SAa '))})
        desc.update({'selection_coeff_homo':int(extract_msms_parameter(first_line, '-Saa '))})

        desc.update({'model':str(self.model_name)})

        # Get the UNIX Time Stamp of when the file was modification
        desc.update({'modification_stamp':os.stat(file_name).st_mtime})

        # Allow deleted files to be tracked in json folder
        desc.update({'active':'active'})

        return desc

    def read_simulations(self, parameter_name='selection_coeff_hetero', max_nrepl=None, verbose=0):
        """
        Read simulations and store into compressed numpy arrays

        Keyword Arguments:
            parameter_name: name of parameter to estimate
            max_nrepl: max nr of replicates per simulated msms file
            verbose: 

        Returns:
            an object of class Genes
        """

        data = []
        positions = []
        description = []

        # Open the directory in which simulation files are stored
        for file_name in os.listdir(self.simulations_folder):

            full_name = self.simulations_folder + '/%s' %(file_name)

            if verbose > 0:
                print(full_name, ': ', end='')

            # Read lines including the metadata
            f = gzip.open(full_name, 'rb')
            file_content = f.read().decode('utf8').split('\n')

            # Search the // char inside the file
            starts = ([i for i, e in enumerate(file_content) if e == '//'])

            # limit the scan to the first max_nrepl items (if set)
            if max_nrepl!=None:
                starts = starts[:max_nrepl]

            if verbose > 0:
                print(len(starts))

            # Populate object with data for each simulated gene
            for idx, pointer in enumerate(starts):

                # Description for each simulation
                description.append(self.extract_description(full_name, file_content[0]))

                nr_columns = int(file_content[pointer+1].split('segsites: ')[1])
                haplotypes = np.zeros((self.nr_samples, nr_columns, 1), dtype='uint8')
                pos = file_content[pointer+2].split(' ')
                pos.pop()
                pos.pop(0)
                positions.append(np.asarray(pos, dtype='float32'))
                del pos

                for j in range(self.nr_samples):

                    hap = list(file_content[pointer + 3 + j])

                    # string processing: if not 0/1 --> convert to 1
                    hap = ['1' if element!='0' and element!=1 else element for element in hap]
                    # switch colours, 1s are black and 0s are white
                    hap = ['255' if element=='1' else element for element in hap]
                    haplotypes[j,:,0] = hap

                data.append(haplotypes)

            f.close()

        gene = ImaGene(data=data, positions=positions, description=description, parameter_name=parameter_name)

        return gene

    def read_VCF(self, verbose=0):
        """
        Read VCF file and store into compressed numpy arrays

        Keyword Arguments:
            verbose: 

        Returns:
            an object of class Genes
        """

        with open(self.VCF_file_name, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]

        header = lines.pop(0)
        ind_pos = header.split('\t').index('POS')
        ind_format = header.split('\t').index('FORMAT')

        nr_individuals = len(header.split('\t')) - ind_format - 1
        nr_sites = len(lines)

        if verbose == 1 | self.nr_samples!=(nr_individuals*2):
            print('Found' + str(nr_individuals) + 'individuals and' + str(nr_sites) + 'sites.')

        haplotypes = np.zeros(((nr_individuals * 2), nr_sites, 1), dtype='uint8')

        data = []
        positions = []
        pos = np.zeros((nr_sites), dtype='int32')

        for j in range(nr_sites):
            # populate genomic position
            pos[j] = int(lines[j].split('\t')[ind_pos])
            # extract genotypes
            genotypes = lines[j].split('\t')[(ind_format+1):]
            genotypes[len(genotypes) - 1] = genotypes[len(genotypes) - 1].split('\n')[0]
            for i in range(len(genotypes)):
                if i == 0:
                    i1 = 0
                    i2 = 1
                else:
                    i2 = i*2
                    i1 = i2 - 1
                if genotypes[i].split('|')[0] == '1':
                    haplotypes[i1,j] = '255'
                if genotypes[i].split('|')[1] == '1':
                    haplotypes[i2,j] = '255'

        positions.append(pos)
        data.append(haplotypes)
        del pos
        del haplotypes

        gene = ImaGene(data=data, positions=positions)

        return gene


class ImaGene:
    """
    A batch of genomic images
    """
    def __init__(self, data, positions, description=[], targets=[], parameter_name=None, classes=[]):
        self.data = data
        self.positions = positions
        self.description = description
        self.dimensions = (np.zeros(len(self.data)), np.zeros(len(self.data)))
        # initialise dimensions to the first image (in case we have only one)
        self.dimensions[0][0] = self.data[0].shape[0]
        self.dimensions[1][0] = self.data[0].shape[1]
        # if reads from real data, then stop here otherwise fill in all info on simulations
        if parameter_name != None:
            self.parameter_name = parameter_name # this is passed by ImaFile.read_simulations()
            self.targets = np.zeros(len(self.data), dtype='int32')
            for i in range(len(self.data)):
                # set targets from file description
                self.targets[i] = self.description[i][self.parameter_name]
                # assign dimensions
                self.dimensions[0][i] = self.data[i].shape[0]
                self.dimensions[1][i] = self.data[i].shape[1]
            self.classes = np.unique(self.targets)
        return None

    def summary(self):
        """
        Prints general info on the object.

        Keyword Arguments:


        Returns:
            0
        """
        nrows = self.dimensions[0]
        ncols = self.dimensions[1]
        print('An object of %d image(s)' % len(self.data))
        print('Rows: min %d, max %d, mean %f, std %f' % (nrows.min(), nrows.max(), nrows.mean(), nrows.std()))
        print('Columns: min %d, max %d, mean %f, std %f' % (ncols.min(), ncols.max(), ncols.mean(), ncols.std()))
        return 0

    def plot(self, index=0):
        """
        Plot one image in gray scale.

        Keyword arguments:
            index: index of image to plot

        Returns:
            0
        """
        image = plt.imshow(self.data[index][:,:,0], cmap='gray')
        plt.show(image)
        return 0

    def majorminor(self):
        """
        Convert to major/minor polarisation.

        Keyword Arguments:

        Returns:
            0
        """
        for i in range(len(self.data)):
            idx = np.where(np.mean(self.data[i][:,:,0]/255., axis=0) > 0.5)[0]
            self.data[i][:,idx,0] = 255 - self.data[i][:,idx,0]
        return 0

    def filter_freq(self, minimal_maf, verbose=0):
        """
        Remove sites whose minor allele frequency is below the set threshold.

        Keyword Arguments:
            minimal_maf: minimal minor allele frequency to retain the site

        Returns:
            0
        """
        for i in range(len(self.data)):
            idx = np.where(np.mean(self.data[i][:,:,0]/255., axis=0) >= minimal_maf)[0]
            self.positions[i] = self.positions[i][idx]
            self.data[i] = self.data[i][:,idx,:]
            # update nr of columns in dimensions
            self.dimensions[1][i] = self.data[i].shape[1]
        return 0

    def resize(self, dimensions=(128, 128), option=None, set_to_boundaries=True):
        """
        Resize all images to same dimensions.

        Keyword Arguments:
            dimensions: tuple, nr of rows and nr of columns
            option: either 'mean', 'min' or 'max'
            set_to_boundaries: if True, all cells are pushed up/down to 255/0

        Returns:
            0
        """
        if option == 'mean':
            dimensions = (int(self.dimensions[0].mean()), int(self.dimensions[1].mean()))
        elif option == 'min':
            dimensions = (int(self.dimensions[0].min()), int(self.dimensions[1].min()))
        elif option == 'max':
            dimensions = (int(self.dimensions[0].max()), int(self.dimensions[1].max()))
        else: pass
        for i in range(len(self.data)):
            image = np.copy(self.data[i][:,:,0])
            self.data[i] = np.zeros((dimensions[0], dimensions[1], 1), dtype='uint8')
            self.data[i][:,:,0] = (skimage.transform.resize(image, dimensions, anti_aliasing=True, mode='reflect')*255).astype('uint8')
            del image
            # reassign data dimensions
            self.dimensions[0][i] = self.data[i].shape[0]
            self.dimensions[1][i] = self.data[i].shape[1]
            if set_to_boundaries == True:
                self.data[i] = (np.where(self.data[i] < 128, 0, 255)).astype('uint8')
        return 0

    def sort(self, ordering):
        """
        Sort rows and/or columns given an ordering.

        Keyword Arguments:
            ordering: either 'rows_freq', 'cols_freq', 'rows_dist', 'cols_dist'

        Returns:
            0
        """
        if ordering == 'rows_freq':
            for i in range(len(self.data)):
                uniques, counts = np.unique(self.data[i], return_counts=True, axis=0)
                counter = 0
                for j in counts.argsort()[::-1]:
                    for z in range(counts[j]):
                        self.data[i][counter,:,:] = uniques[j,:,:]
                        counter += 1
        elif ordering == 'cols_freq':
            for i in range(len(self.data)):
                uniques, counts = np.unique(self.data[i], return_counts=True, axis=1)
                counter = 0 #
                for j in counts.argsort()[::-1]:
                    for z in range(counts[j]):
                        self.data[i][:,counter,:] = uniques[:,j,:]
                        counter += 1
        elif ordering == 'rows_dist':
            for i in range(len(self.data)):
                uniques, counts = np.unique(self.data[i], return_counts=True, axis=0)
                # most frequent row in float
                top = uniques[counts.argsort()[::-1][0]].transpose().astype('float32')
                # distances from most frequent row
                distances = np.mean(np.abs(uniques[:,:,0] - top), axis=1)
                # fill in from top to bottom
                counter = 0
                for j in distances.argsort():
                    for z in range(counts[j]):
                        self.data[i][counter,:,:] = uniques[j,:,:]
                        counter += 1
        elif ordering == 'cols_dist':
            for i in range(len(self.data)):
                uniques, counts = np.unique(self.data[i], return_counts=True, axis=1)
                # most frequent column
                top = uniques[:,counts.argsort()[::-1][0]].astype('float32')
                # distances from most frequent column
                distances = np.mean(np.abs(uniques[:,:,0] - top), axis=0)
                # fill in from left to right
                counter = 0
                for j in distances.argsort():
                    for z in range(counts[j]):
                        self.data[i][:,counter,:] = uniques[:,j,:]
                        counter += 1
        else:
            print('Select a valid ordering.')
            return 1
        return 0

    def convert(self, normalise=False, flip=False, verbose=False):
        """
        Check for correct data type and convert otherwise. Convert to float numpy arrays [0,1] too. If flip true, then flips 0-1
        """
        # if list, put is as numpy array
        if type(self.data) == list:
            if len(np.unique(self.dimensions[0]))*len(np.unique(self.dimensions[1])) == 1:
                if verbose:
                    print('Converting to numpy array.')
                self.data = np.asarray(self.data)
            else:
                print('Aborted. All images must have the same shape.')
                return 1
        # if unit8, put it as float and divide by 255
        if self.data.dtype == 'uint8':
            if verbose:
                print('Converting to float32.')
            self.data = self.data.astype('float32')
        if self.data.max() > 1:
            if verbose:
                print('Converting to [0,1].')
            self.data /= 255.
        # normalise
        if normalise==True:
            if verbose:
                print('Normalising samplewise.')
            for i in range(len(self.data)):
                mean = self.data[i].mean()
                std = self.data[i].std()
                self.data[i] -= mean
                self.data[i] /= std
        # flip
        if flip==True:
            if verbose:
                print('Flipping values.')
            for i in range(len(self.data)):
                self.data[i] = 1. - self.data[i]
        if verbose:
            if self.data.shape[0] > 1: 
                print('A numpy array with dimensions', self.data.shape, 'and', len(self.targets), 'targets and', len(self.classes), 'classes.')
            else: # one real image
                print('A numpy array with dimensions', self.data.shape)
        return 0

    def set_classes(self, classes=[], nr_classes=0):
        """
        Set classes (or reinitiate)
        """
        # at each call reinitialise for safety
        targets = np.zeros(len(self.data), dtype='int32')
        for i in range(len(self.data)):
            # set target from file description
            targets[i] = self.description[i][self.parameter_name]
        self.classes = np.unique(targets)
        # calculate and/or assign new classes
        if nr_classes > 0:
            self.classes = np.asarray(np.linspace(targets.min(), targets.max(), nr_classes), dtype='int32')
        elif len(classes)>0:
            self.classes = classes
        del targets
        return 0

    def set_targets(self):
        """
        Set targets for binary or categorical classification (not for regression) AFTER running set_classes
        """
        # initialise
        self.targets = np.zeros(len(self.data), dtype='int32')
        for i in range(len(self.targets)):
            # reinitialise
            self.targets[i] = self.description[i][self.parameter_name]
            # assign label as closest class
            self.targets[i] = self.classes[np.argsort(np.abs(self.targets[i] - self.classes))[0]]
        return 0

    def subset(self, index):
        """
        Subset object to index array (for shuffling or only for multiclassification after setting classes and targets)
        """
        # update based on index
        self.targets = self.targets[index]
        self.data = self.data[index]
        self.positions = [self.positions[i] for i in index]
        self.description = [self.description[i] for i in index]
        for i in range(len(self.data)):
            self.dimensions[0][i] = self.data[i].shape[0]
            self.dimensions[1][i] = self.data[i].shape[1]
        return 0

    def save(self, file):
        """
        Save to file
        """
        with open(file, 'wb') as fp:
            pickle.dump(self, fp)
        return 0

    def crop(self, window):
        """
        crop or extend haplotype window for genomic image object. Window size are adjusted from center

        Arguments:
            window: haplotype window size

        """

        for i, image in enumerate(self.data):
            x, y, c = image.shape[0], image.shape[1], image.shape[2]

            if y == window:
                continue
            
            #when even no. haplotype column
            if y % 2 == 0:
                if window < y:
                    starty = y // 2 - window // 2
                    self.data[i] = image[:, starty:starty + window, :]

                #perform padding
                else:
                    padding_len = (window - y) // 2
                    padding = np.zeros((x, padding_len, c))
                    self.data[i] = np.concatenate((padding, image, padding), axis=1)
            
            #when odd no.haplotype column
            #will result in slight offset for window by padding a empty padding on the right hand side
            else:
                offset_padding = np.zeros((x, 1, c))
                image = np.concatenate((image, offset_padding), axis = 1)
                #perform cropping
                if window < y:
                    starty = y // 2 - window // 2
                    self.data[i] = image[:, starty:starty + window, :]

                #perform padding
                else:
                    padding_len = (window - y) // 2
                    padding = np.zeros((x, padding_len, c))
                    self.data[i] = np.concatenate((padding, image, padding), axis=1)


            #update dimension
            self.dimensions[0][i] = self.data[i].shape[0]
            self.dimensions[1][i] = self.data[i].shape[1]

        return None

            
                


class ImaNet:
    """
    Training and Learning
    """
    def __init__(self, name=None, model=None):
        self.name = name
        self.scores = {'val_loss': [], 'val_accuracy': [], 'loss': [], 'accuracy': [], 'mae': [], 'val_mae': []}
        self.test = np.zeros(2)
        self.values = None # matrix(3,nr_test) true, map, mle
        return None

    def update_scores(self, score):
        """
        Append new scores after each training
        """
        for key in self.scores.keys():
            if key in score.history:
                self.scores[key].append(score.history[key])
        return 0

    def plot_train(self, file=None):
        """
        Plot training accuracy/mae and loss/mse
        """
        loss = self.scores['loss']
        val_loss = self.scores['val_loss']
        # if regression
        if len(self.scores['mae'])>0:
            acc = self.scores['mae']
            val_acc = self.scores['val_mae']
            label = 'mae'
        else: # if not
            acc = self.scores['accuracy']
            val_acc = self.scores['val_accuracy']
            label = 'accuracy'
        epochs = range(1, len(loss) + 1)

        plt.figure()
        plt.subplots_adjust(wspace = 0, hspace = 0.4)
        plt.subplot(211)

        plt.plot(epochs, loss, 'bo', label='Training loss')
        plt.plot(epochs, val_loss, 'b', label='Validation loss')
        plt.title('Training and validation loss')
        plt.legend()

        plt.subplot(212)

        plt.plot(epochs, acc, 'bo', label='Training '+label)
        plt.plot(epochs, val_acc, 'b', label='Validation '+label)
        plt.title('Training and validation '+label)
        plt.legend()

        if file==None:
            plt.show()
        else:
            plt.savefig(file)

        return 0

    def predict(self, gene, model):
        """
        Calculate predicted values (many, I assume this is for testing not for single prediction); output is a matrix with rnows=2, row 0 is true, row 1 is MAP, row 2 is posterior mean
        """
        self.values = np.zeros((3, gene.data.shape[0]), dtype='float32')
        # if binary or regression
        if len(gene.targets.shape) == 1:
            probs = model.predict(gene.data, batch_size=None)[:,0]
            self.values[1,:] = np.where(probs < 0.5, 0., 1.)
            self.values[0,:] = gene.targets
            self.values[2,:] = probs
        else:
            probs = model.predict(gene.data, batch_size=None)
            self.values[1,:] = gene.classes[np.argmax(probs, axis=1)]
            self.values[0,:] = gene.classes[np.argmax(gene.targets, axis=1)]
            self.values[2,:] = [np.average(gene.classes, weights=probs[i]) for i in range(probs.shape[0])]

        return 0

    def plot_scatter(self, MAP=True, file=None):
        """
        Plot scatter plot (on testing set)
        """
        # if MAP
        if MAP == True:
            plt.scatter(self.values[0,:], self.values[1,:], marker='o')
        else: # if regression
            plt.scatter(self.values[0,:], self.values[2,:], marker='o')
        #plt.title('Relationship between true and predicted values')
        plt.xlabel('True')
        plt.ylabel('Predicted')
        if file==None:
            plt.show()
        else:
            plt.savefig(file)
            plt.close()
        return 0

    def plot_cm(self, classes, file=None, text=False):
        """
        Plot confusion matrix (on testing set)
        """
        cm = confusion_matrix(self.values[0,:], self.values[1,:])
        accuracy = np.trace(cm) / float(np.sum(cm))
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

        fig = plt.figure()
        plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
        plt.title('Confusion matrix')
        plt.colorbar()

        tick_marks = np.arange(len(classes))
        plt.xticks(tick_marks, classes, fontsize=8)
        plt.yticks(tick_marks, classes, fontsize=8)

        thresh = cm.max() / 1.5
        if text==True:
            for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
                plt.text(j, i, "{:0.4f}".format(cm[i, j]), horizontalalignment="center", color="white" if cm[i, j] > thresh else "black")

        plt.ylabel('True label')
        plt.xlabel('Predicted label\naccuracy={:0.4f}'.format(accuracy))
        plt.tight_layout()

        if (file==None):
            plt.show()
        else:
            plt.savefig(file)
            plt.close()

        return 0

    def save(self, file):
        """
        Save to file
        """
        with open(file,'wb') as fp:
            pickle.dump(self, fp)

        return 0





