# GraSSRep

This repository contains all the necessary code and commands for "GraSSRep: Graph-Based Self-Supervised Learning for Repeat Detection in Metagenomic Assembly."

We propose GraSSRep, a novel approach that leverages the assembly graph’s structure through graph neural networks (GNNs) within a self-supervised learning framework to classify DNA sequences (contigs) into repetitive and non-repetitive categories. Here is an overview of our method:

![Overview](https://github.com/aliaaz99/GraSSRep/assets/136205616/43a06c8b-9f53-4359-b47d-b7abf4c5e36b)


## Installation & Dependencies

The code is based on Python 3.7 and should run on Unix-like operating systems (MacOS, Linux).

To install the dependencies for GraSSRep, you can use the `environment.yml` file provided to build a Conda environment with all necessary dependencies. 
The GraSSRep environment will need to be activated for each usage. 

```sh
conda env create -f environment.yml
conda activate GraSSRep
```

Additionally, you will need to install other packages using `pip` after creating and activating the environment.

```sh
pip install -r requirements. txt
```

Finally, make sure you have `Spades` and `MUMmer` installed. You can follow the installation instructions provided in the [Spades](https://github.com/ablab/spades) and [MUMmer](https://mummer.sourceforge.net/) repository.



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

For example, for the shakya_1 dataset:

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

For simulated data, run `SimulatedData.sh` to generate reference genomes, read, assemble them, and run the code for the three different cases discussed in the paper.

2. **Generating the results:**

For any dataset with the read pairs and possible reference genome, there are two main steps to detect the repeated contigs:

  i) **Assembly Graph:**
  
In this step, reads are assembled into contigs, the assembly graph is provided by `spades`, and the features are calculated for the contigs or nodes of the assembly graph. Note that if the reference genome is available, reads are generated from the reference genome using `wgsim `. Ground truth labels for repeats and non-repeat contigs are also found based on the reference genome using `Nucmer`.
All these steps are executed in the `mainGraph.py` code. In other words, this code includes Steps 1 and 2 in the overview figure. 
You can vary the following parameters:

| Parameter          | Default Value     | Description                                                              |
|--------------------|-------------------|--------------------------------------------------------------------------|
| `--name`           | "shakya_1"     | Folder containing reads and reference genomes (if available) in Data folder |
| `--read1`          | "outRead1.fq"     | Read 1 file name                                                         |
| `--read2`          | "outRead2.fq"     | Read 2 file name                                                         |
| `--assembly`       | 1                 | Flag to run the metaSpades or not                                              |
|                    |                | (if the assembly graph or spades is generated before, change this to 0)           |
| `--idy`            | 95                | Identity for repeat detection in %                                       |
| `--nL`             | 95                | Normalized length for repeat detection in %                              |
| `--cN`             | 2                 | Copy number for repeat detection                                         |
| `--isGT`           | 1                 | Flaf for the availability of the ground truth (reference genome)                      |
|                    |                   | (if there is no reference genome available, change this to 0)            |
| `--num_processes`  | 30                | Number of processors                                                     |


Please note that the reference genome is used solely for finding the ground truth contigs and evaluating the model. If a dataset lacks a reference genome, the code can be modified to utilize the provided read pairs for identifying repetitive contigs.

  ii) **Repeat detection:**

Finally, using the assembly graph structure and node features, the `mainRepeatDetection.py` code classifies contigs into repeat and non-repeat classes using GNN and in a self-supervised manner and finally applies a fine tuning step. This includes Steps 3, 4, and 5 in the overview figure. Parameters that you can vary for this code are:

| Parameter            | Default Value | Description                                                              |
|----------------------|---------------|--------------------------------------------------------------------------|
| `--name`             | "shakya_1"    | Folder containing reads and reference genomes (if available) in Data folder |
| `--p`                | 35            | p threshold to generate pseudo-labels and apply fine-tuning              |
| `--N_iter`           | 10            | Number of iterations to repeat the results and average it                |
| `--isSemiSupervised`| 0             | Flag to use semi-supervised learning instead of self-supervised learning |
| `--noGNN`            | 0             | Flag to exclude GNN step or not                                          |
| `--noRF`             | 0             | Flag to exclude RF step or not                                           |


This three-step process reads data files from the `Data` directory, generates results for each setup specified in the script, and saves them in the corresponding folder within the `Results` directory.

The output file `final_pred.csv` contains the unitigs' names, their indices in the unitig graph, along with true labels and predicted labels.




