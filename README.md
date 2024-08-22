# ImaGene

**ImaGene** is a supervised machine learning algorithm to predict natural selection and estimate selection coefficients from population genomic data.
It can be used to estimate any parameter of interest from an evolutionary population genetics model.

**ImaGene** implements a convolutional neural network (CNN) which takes as input haplotypes of a _locus_ of interest for a population.
It outputs confusion matrices as well as point estimates of the selection coefficient (or any parameter of interest) along with its posterior distribution and various metrics of confidence.

### Citation

The original manuscript can be found [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2927-x) and it is open access.
You should cite it as:

Torada, L., Lorenzon, L., Beddis, A. _et al_. ImaGene: a convolutional neural network to quantify natural selection from genomic data. _BMC Bioinformatics_ __20__, 337 (2019)

doi:10.1186/s12859-019-2927-x

and you can download the citation file in [.ris](citeme.ris) or [.json](citeme.json) format.

Related studies from the group are on [balancing selection](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13379), [adaptive introgression](https://elifesciences.org/articles/64669), and a [review](https://doi.org/10.1093/gbe/evad008) paper.

### Download and installation

We provide examples on how to run **ImaGene** on your local machine, on Google Colab or Google Cloud in the `Tutorials` folder.
Below are instructions to download and install **ImaGene** on your local machine or server.

Download the repository using git.
```
git clone https://github.com/mfumagalli/ImaGene
```

**ImaGene** runs under Python3 and it is interfaced with [tensorflow](https://www.tensorflow.org) and [keras](https://keras.io/).
We recommend using [mamba](https://mamba.readthedocs.io/en/latest/) to set the environment and take care of all dependencies.
There are detailed instructions on how to download conda for [linux](https://conda.io/docs/user-guide/install/linux.html) and [macOS](https://conda.io/docs/user-guide/install/macos.html).
A suitable environment can be created with

`mamba create -n ImaGene python=3 tensorflow=2 keras=2 numpy scipy scikit-image scikit-learn matplotlib pydot arviz jupyterlab -y`

Activate the environment with 

`conda activate ImaGene`

and pip-install _protobuf_ with

`pip3 install --force --no-deps protobuf`

You can deactivate the environment with

`conda deactivate`.


**ImaGene** receives training data in _msms_.
Unless you have already generated training data, you are required to download _msms_ separately following the instructions [here](https://www.mabs.at/publications/software-msms/downloads/).
Follow the link, download the .zip folder and extract it.
The .jar file of interest will be in the `lib` folder.
There are no requirements for msms to be installed in a specific folder.
However, _msms_ requires java to be installed.
On unix Debian systems just type `sudo apt-get update && apt-get upgrade; sudo apt-get install default-jdk`
Otherwise follow the link [here](https://www.java.com/en/download/) if you need to install java on other systems.
Remember that java must be in your /usr/bin folder.
In unix systems you can create a symbolic link with `ln -s ~/Downloads/java-XXX/jre/bin/java /usr/bin/java`, as an example.

### Usage

Please look at the jupyter notebook `01_binary.ipynb` (or the corresponding `Colab` version) for a tutorial on how to use **ImaGene** for predicting natural selection with a simple binary classification.
We also provide examples on how **ImaGene** can be used for multiclass classification in `02_multiclass.ipynb` and `03_multiclass_for_continuous.ipynb` and for regression in `04_regression.ipynb` (or the corresponding `Colab` versions).

Finally, we provide an utility `generate_dataset.sh` to quickly generate simulations with _msms_ to be used for training. 
This script takes an input file with all parameters needed for the simulations.
An example of this file is `params.txt` and tutorials show how to run such simulations in practice.
More information can be found in the tutorials.

The folder `Reproduce` contains all scripts used for the analyses shown in the manuscript.

### Contributors (in alphabetical order)

- main: [Matteo Fumagalli](https://www.qmul.ac.uk/sbbs/staff/matteo-fumagalli.html)
- others (in alphabetical order): Alice Beddis, Ulas Isildak, Sirimar (Nook) Laosinwattana, Lucrezia Lorenzon, Jacky Pui Chung Siu, Luis Torada


