# ImaGene

**ImaGene** implements a supervised machine learning algorithm to predict natural selection and estimate selection coefficients from population genomic data.
Specifically, it uses a convolutional neural network (CNN) which takes as input haplotypes for a population and locus of interest.
It outputs confusion matrcies as well as point estimates of the selection coefficient along with its posterior distribution and various metrics of confidence.

Download the repository using git.
```
git clone https://github.com/mfumagalli/ImaGene
```

**ImaGene** runs under Python3 and it inferfaces with keras.
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

Please look at the jupyter notebook "Quick start" for a short tutorial on using **ImaGene** for predicting natural selection.

The folder Reproduce contains all scripts used for the analyses shown in the manuscript.

### Contributors (in alphabetical order)
Alice Beddis, Matteo Fumagalli, Ulas Isildak, Lucrezia Lorenzon, Luis Torada



