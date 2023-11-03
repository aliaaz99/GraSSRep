# GraSSRep

This repository contains all the necessary code and commands for "GraSSRep: Graph-Based Self-Supervised Learning for Repeat Detection in Metagenomic Assembly."

We propose GraSSRep, a novel approach that leverages the assembly graph’s structure through graph neural networks (GNNs) within a self-supervised learning framework to classify DNA sequences into repetitive and non-repetitive categories


## Installation

The code is based on Python 3.7 and should run on Unix-like operating systems (MacOS, Linux).

### Python libraries

Make sure you have the Python packages listed in `requirements.txt` installed. You can install them using the following command:

```sh
$ pip install -r requirements.txt
```

### Packages:

In addition, ensure that you have installed these required packages:

- **Wgsim**: Follow the installation instructions provided in the [wgsim repository](https://github.com/lh3/wgsim).
- **ABySS**: Follow the installation instructions provided in the [Abyss repository](https://github.com/bcgsc/abyss).
- **Bowtie2**: Follow the installation instructions provided on the [bowtie2 website](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
- **{MUM}mer**: Follow the installation instructions provided in the [mummer repository](https://mummer4.github.io/).


## Running the Codes

1. **Folder Structure:**

You should have the following directory structure in the project folder:

```
├── Codes.py
├── Results
├── Data
│ ├── shakya_1
│ │ ├── (read pairs and reference genome)
│ ├── shakya_2
│ │ ├── (read pairs and reference genome)
│ ├── (other folders for different cases)
│ │ ├── (read pairs and reference genome)
```

You need to place your data files, including read pairs, in `.fq` format and reference genome in `.fasta` format in the respective folders inside the `Data` directory.

For example for the shakya_1 dataset:

```
├── Data
│ ├── ...
│ ├── shakya_2
│ ├── shakya_1
│ │ ├── outRead1.fq
│ │ ├── outRead2.fq
│ │ ├── ref_genome.fasta
│ ├── ...
```

For simulated data, run `genSimData.py` to generate reference genomes for the three different cases discussed in the paper. Read pairs will be generated in the next steps.

2. **Generating the results:**

For any dataset with the read pairs and a reference genome, there are three main steps:

  i) **Assembly**
In this step, reads are generated from the reference genome (if available), assembled into unitigs, and mapped to the unitigs. Ground truth for repeats and non-repeat unitigs is calculated based on the reference genome.
All these steps are executed in the `sequencing.py` code.
You can vary the number of reads to control coverage.

Please note that the reference genome is used solely for calculating the ground truth and evaluating the model. If a dataset lacks a reference genome, the code can be modified to utilize the provided read pairs for identifying repetitive unitigs.

  ii) **Unitig graph construction and feature extraction**
Using read mapping information, the `unitigGraphFeatures.py` code generates and saves the unitig graph with node observations for each node on the graph. Criteria for ground truth repeats can be adjusted in this code.

  iii) **Repeat detection**
Finally, using the unitig graph and node features, the `repeatDetection.py` code classifies unitigs into repeat and non-repeat classes. The threshold $p$ can be varied in this code.

This three-step process reads data files from the `Data` directory, generates results for each setup specified in the script, and saves them in the corresponding folder within the `Results` directory.

The output file `final_pred.csv` contains the unitigs' names, their indices in the unitig graph, along with true labels and predicted labels.




