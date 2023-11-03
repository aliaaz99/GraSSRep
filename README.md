# Graph-Based Self-Supervised Learning for Repeat Detection in Metagenomic Assembly

This repository contains all necessary code and commands to detect the repeated sequences in a genomic dataset.

Here we propose a novel approach that leverages the assembly graph’s structure through graph neural networks (GNNs) within a self-supervised
learning framework to classify DNA sequences (unitigs) into repetitive and non-repetitive categories.
We frame this problem as a node classification task within the assembly graph.


## Installation

The code is based on Python 3.7 and should run on Unix-like operating systems (MacOS, Linux).

### Python libraries

Make sure you have the python packages listed in `requirements.txt` installed. You can install them using the following command:

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
├──Codes.py
├──Results
├── Data
│ ├── simulated_L400_C25
│ │ ├── (read pairs and reference genome)
│ ├── shakya_1
│ │ ├── (read pairs and reference genome)
│ ├── (other folders for different cases)
│ │ ├── (read pairs and reference genome)
```

You need to place your data files, including read pairs, in `.fq` format and reference genome in `.fasta` format in the respective folders inside the `Data` directory.

You need to have three main files provided, for example for shakya_1 dataset:

```
├── Data
│ ├── shakya_1
│ │ ├── outRead1.fq
│ │ ├── outRead2.fq
│ │ ├── ref_genome.fasta
```

2. **Getting the results:**

For any dataset with the read pairs and a reference genome, there are three main steps:

  1. **Assembly**
In this step, reads are generated from the reference genome (if we have the data), then the reads are assembled into the unitigs. After that, reads are mapped to the unitigs, and the mapping data is saved. Also, the ground truth for the repeats and non-repeat unitigs is calculated based on the reference genome. All these steps are done in the `sequencing.py` code.

Note that the reference genome is just used to calculate the ground truth and evaluate the model. So if in a dataset there is no reference genome, the code can be edited to just use the provided read pairs and find the repetitive unitgis. 

  3. **Unitig graph construction and feature extraction**
Using the read mapping information from the previous part, the `unitigGraphFeatures.py` code generates and saves the unitig graph with the node observations for each node on the graph.

  4. **Repeat detection**
Finally, using the unitig graph and the features for each node, `repeatDetection.py` code classifies the unitigs to repeat and non-repeat classes.



This three step process will read the data files located in the `Data` directory, generate results for each setup specified in the script, and save them in the corresponding folder in the `Results` folder.

The output file `final_pred.csv` have the final results. It has the name of the unitigs, with their index in the unitig graph, along with the true labels and predicted labels.




