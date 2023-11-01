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
All necessary dependencies are listed in requirements.txt. You can install them with:

```sh
$ pip install -r requirements.txt
```

In addition, you will need wgsim, ABySS, Bowtie2, and Mummer.
