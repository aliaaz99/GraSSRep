# GraSSRep

This repository contains all the necessary code and commands for "GraSSRep: Graph-Based Self-Supervised Learning for Repeat Detection in Metagenomic Assembly."

We propose GraSSRep, a novel approach that leverages the assembly graph’s structure through graph neural networks (GNNs) within a self-supervised learning framework to classify DNA sequences (contigs) into repetitive and non-repetitive categories. Here is an overview of our method:

![Overview](https://github.com/aliaaz99/GraSSRep/assets/136205616/43a06c8b-9f53-4359-b47d-b7abf4c5e36b)


## Installation & Dependencies

The code is based on Python 3.7 and should run on Unix-like operating systems (MacOS, Linux).

Dependencies for GraSSRep include: 

You can use the `environment.yml` file provided to build a Conda environment with all necessary dependencies. 
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

For any dataset with the read pairs and a reference genome, there are two main steps:

  i) **Assembly:**
  
In this step, reads are assembled into contigs, and then the assembly graph is provided by `spades`. Note that if the reference genome is available, reads are generated from the reference genome using `wgsim `. Ground truth labels for repeats and non-repeat contigs are also found based on the reference genome using `Nucmer`.
Also, 
All these steps are executed in the `mainGraph.py` code.
You can vary the following parameters:

| Parameter          | Description                                                              | Default Value     |
|--------------------|--------------------------------------------------------------------------|-------------------|
| `--name`           | Folder containing reads and reference genomes (if available) in Data folder |                   |
| `--read1`          | Read 1 file name                                                         | "outRead1.fq"    |
| `--read2`          | Read 2 file name                                                         | "outRead2.fq"    |
| `--assembly`       | Apply the metaSpades or not                                              | 1                 |
|                     (if the assembly graph is generated before, change this to 0)                              |
| `--idy`            | Identity for repeat detection in %                                       | 95                |
| `--nL`             | Normalized length for repeat detection in %                              | 95                |
| `--cN`             | Copy number for repeat detection                                         | 2                 |
| `--isGT`           | Availability of the ground truth (reference genome)                      | 1                 |
|                    | (if there is no reference genome available, change this to 0)            |                   |
| `--num_processes`  | Number of processors                                                     | 30                |



Please note that the reference genome is used solely for finding the ground truth contigs and evaluating the model. If a dataset lacks a reference genome, the code can be modified to utilize the provided read pairs for identifying repetitive contigs.

  ii) **Repeat detection:**

Finally, using the unitig graph and node features, the `repeatDetection.py` code classifies unitigs into repeat and non-repeat classes. The threshold $p$ can be varied in this code.

This three-step process reads data files from the `Data` directory, generates results for each setup specified in the script, and saves them in the corresponding folder within the `Results` directory.

The output file `final_pred.csv` contains the unitigs' names, their indices in the unitig graph, along with true labels and predicted labels.




