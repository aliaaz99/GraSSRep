# Graph-Based Self-Supervised Learning for Repeat Detection in Metagenomic Assembly

This repository contains all necessary code and commands to detect the repeated sequences in a genomic dataset.

Here we propose a novel approach that leverages the assembly graph’s structure through graph neural networks (GNNs) within a self-supervised
learning framework to classify DNA sequences (unitigs) into repetitive and non-repetitive categories.
We frame this problem as a node classification task within the assembly graph.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)


## Installation

The code is based on Python 3.7 and should run on Unix-like operating systems (MacOS, Linux).

### Python libraries

Make sure you have the following dependencies installed:

- Python packages listed in `requirements.txt` file. You can install them using the following command:

```sh
$ pip install -r requirements.txt
```

### Packages:

In addition, ensure that you have installed these required packages:

- [wgsim](https://github.com/lh3/wgsim): Follow the installation instructions provided in the [wgsim repository](https://github.com/lh3/wgsim).
- [Abyss](https://github.com/bcgsc/abyss): Follow the installation instructions provided in the [Abyss repository](https://github.com/bcgsc/abyss).
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): Follow the installation instructions provided on the [bowtie2 website](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
- [mummer](https://mummer4.github.io/): Follow the installation instructions provided in the [mummer repository](https://mummer4.github.io/).



