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
pip install -r requirements.txt
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

For simulated data, you can run `SimulatedData.sh` to simulate reference genomes, generate read-pairs, assemble them, and run the code for the three different cases discussed in the paper.
```sh
chmod +x SimulatedData.sh
./SimulatedData.sh
```

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
|                    |                | (if the assembly graph is generated before, change this to 0)           |
| `--idy`            | 95                | Identity for ground truth repeat detection in %                                       |
| `--nL`             | 95                | Normalized length for ground truth repeat detection in %                              |
| `--cN`             | 2                 | Copy number for ground truth repeat detection                                         |
| `--isGT`           | 1                 | Flaf for the availability of the reference genome                     |
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
| `--gnn_dim_hidden`   | "16,16"       | Hidden units per layer for GNN, comma seperated                          |
| `--isSemiSupervised` | 0             | Flag to use semi-supervised learning instead of self-supervised          |
| `--noGNN`            | 0             | Flag to exclude GNN step or not                                          |
| `--noRF`             | 0             | Flag to exclude RF step or not                                           |


## Output files:

The above two-step process reads data files from the `Data` directory, generates results for each setup specified in the script, and saves them in the corresponding folder within the `Results` directory. 
The main files generated by the project are listed below:

- **before_rr.fasta**: Contains contigs (nodes) from the assembly graph provided by SPAdes.
- **assembly_graph.fastg**: Represents the connections between contigs, forming the assembly graph structure.

  *Note:* The above two files, provided by SPAdes, result from read-pairs and serve as the primary inputs for our method. If these files are present in any other dataset, we can exclude assembly with SPAdes from the process. The following files are generated by GraSSRep:

- **Bandage_map.csv**: Provides the mapping between contig labels used for plotting with Bandage and the labels used in our method, in sequential order.
- **G_node_to_index.csv**: Contains three columns:
  - `contig`: Represents the label of contigs in sequential order as mentioned in the previous file.
  - `idx`: Denotes the index of each contig on the graph for Graph Neural Network (GNN).
  - `y_true_B`: Indicates the true label of each contig.
- **final_pred.csv**: Shares the same columns as the previous file, with two additional columns:
  - `y_pred_B`: Represents the predicted label for each contig.
  - `train_mask`: Indicates whether each contig was included in the training set or not.

## Example usage:

Here is an example workflow demonstrating how to use GraSSRep. The following commands generate a reference genome, produce reads, create an assembly graph, and detect repeats:

### Generating reference genome:
In this step, we generate one of the cases of the Simulated dataset introduced in the Results section. The `insertRepeat.py` script is used to generate a reference genome with inserted repeats. The `--name` parameter specifies the name of the dataset, and the `--L` and `--C` parameters specify the length and copy number of the inserted repeat sequences. A directory for the results is created using mkdir to store the results in the following steps.
```sh
python insertRepeat.py --name Example --L 400 --C 25
mkdir -p Results/Example
```
Note that if the reference genomes are available, there is no need to generate them, and you can skip this step by just placing your reference genome as ref_genome.fasta in the desired folder. Also, if there is no reference genome available, you can skip this step and simply continue from the next step.

### Generating the reads:
The wgsim tool simulates read pairs from the generated reference genome in the previous step. 

```sh
wgsim -N 2000000 -1 101 -2 101 "Data/Example/ref_genome.fasta" "Data/Example/outRead1.fq" "Data/Example/outRead2.fq" > "Results/Example/wgsim.log"
```

The parameters -N, -1, and -2 specify the number of reads and the lengths of the paired-end reads, respectively. The outputs are saved as outRead1.fq and outRead2.fq in the specified data directory. The log file for this step is saved in the results directory.

### Generate assembly graph:
The `mainGraph.py` script assembles the reads into contigs and generates the assembly graph as explained previously.

```sh
python mainGraph.py --name Example --read1 outRead1.fq --read2 outRead2.fq --assembly 1 --idy 95 --nL 95 --cN 2 --isGT 1 > "Results/Example/assemblygraph.log"
```
### Detecting repeats:

The `mainRepeatDetection.py` script uses the assembly graph structure and node features to classify contigs into repeat and non-repeat categories. 

```sh
python mainRepeatDetection.py --p 50 --name "Example" --N_iter 5 --gnn_dim_hidden 16,16 > "Results/Example/repeatDetection.log"
```

The `--p` parameter specifies the threshold for generating pseudo-labels and applying fine-tuning. The log file for this step is saved in the results directory.

In addition to the repeat detection results, other files containing information on the contigs are stored as explained previously in the Output Files section.
