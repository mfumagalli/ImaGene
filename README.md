# ImaGene

**ImaGene** implements a supervised machine learning algorithm to predict natural selection and estimate selection coefficients from population genomic data.
Specifically, it uses a convolutional neural network (CNN) which takes as input haplotypes for a population and locus of interest.
It outputs confusion matrices as well as point estimates of the selection coefficient along with its posterior distribution and various metrics of confidence.

### Citation

The original manuscript can be found [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2927-x) and it is open access.
You should cite it as:

Torada, L., Lorenzon, L., Beddis, A. _et al_. ImaGene: a convolutional neural network to quantify natural selection from genomic data. _BMC Bioinformatics_ __20__, 337 (2019)

doi:10.1186/s12859-019-2927-x

and you can download the citation file in [.ris](citeme.ris) or [.json](citeme.json) format.

### Download and installation

Download the repository using git.
```
git clone https://github.com/mfumagalli/ImaGene
```

**ImaGene** runs under Python3 and it is inferfaced with [keras](https://keras.io/).
We recommend using [conda](https://conda.io/docs/index.html) to set the environment and take care of all dependencies.
There are detailed instructions on how to download conda for [linux](https://conda.io/docs/user-guide/install/linux.html) and [macOS](https://conda.io/docs/user-guide/install/macos.html).

**ImaGene** is currently interfaced with [msms](https://www.mabs.at/ewing/msms/index.shtml) but you are required to download it separately following the instructions [here](https://www.mabs.at/ewing/msms/download.shtml).
Follow the link, download the .zip folder and extract it.
The .jar file of interest will be in the `lib` folder.
There are no requirements for msms to be installed in a specific folder.
However, msms requires java to be installed.
On unix Debian systems just type `sudo apt-get update && apt-get upgrade; sudo apt-get install default-jdk`
Otherwise follow the link [here](https://www.java.com/en/download/) if you need to install java on other systems.
Remember that java must be in your /usr/bin folder.
In unix systems you can create a symbolic link with `ln -s ~/Downloads/java-XXX/jre/bin/java /usr/bin/java`, as an example.

### Usage

Please look at the jupyter notebook `binary.ipynb` for a short tutorial on
how to use **ImaGene** for predicting natural selection with a simple binary classification.
We also provide examples on how **ImaGene** can be used for multiclass classification in `multiclass.ipynb` and `continuous.ipynb`.

Finally, we provide an utility `generate_dataset.sh` to quickly generate simulations with msms to be used for training. 
This script takes an input file with all parameters needed for the simulations.
An example of this file is `params.txt` and tutorials show how to run such simulations in practice.

The folder `Reproduce` contains all scripts used for the analyses shown in the manuscript.
The folder `HPC` within `Reproduce` should be ignored.

### Contributors (in alphabetical order)
Alice Beddis, Matteo Fumagalli, Ulas Isildak, Lucrezia Lorenzon, Luis Torada


