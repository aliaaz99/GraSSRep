from toolsGraph import *
import sys
import time
import argparse as ap
import os
import json
 # """# **Initialize**"""

parser = ap.ArgumentParser()
parser.add_argument('--name', help='Folder containing reads and reference genomes (if available) in Data folder', default="shakya_1")
parser.add_argument('--read1', help='Read 1 file name', default="outRead1.fq")
parser.add_argument('--read2', help='Read 2 file name', default="outRead2.fq")
parser.add_argument('--assembly', type=int, help='Apply the metaSpades or not', default=1)
parser.add_argument('--idy', type=int, help='Identity for repeat detection in %', default=95)
parser.add_argument('--nL', type=float, help='Normalized length for repeat detection in %', default=95)
parser.add_argument('--cN', type=int, help='Copy number for repeat detection', default=2)
parser.add_argument('--isGT', type=int, help='Availability of the ground truth (reference genome)', default=1)
parser.add_argument('--num_processes', type=int, help='Number of processors', default=30)
args = parser.parse_args()

# name of files:
read1 = 'Data/' + args.name + '/' + args.read1
read2 = 'Data/' + args.name + '/' + args.read2
graphContigs = 'Data/' + args.name + '/assemblies/assembly_graph.fastg'
before_rr = 'Data/' + args.name + '/assemblies/before_rr.fasta'
coords_file_name = 'Data/' + args.name + '/assemblies/contig_mapping.coords'

# Creat folders to save results of data
try:
    os.mkdir('Results/' + args.name)
except:
    print("Results folder already exists")   


# save the parameters of the run:
with open('Results/' + args.name + '/readmapping_parameters.txt', 'w+') as f:
    json.dump(args.__dict__, f, indent=2)

# Run the pipeline and save the outputs in a file:
with open('Results/' + args.name + '/readmapping_out.txt', "w") as f:
    sys.stdout = f

    start_time = time.time()
    """Assemble the reads and map contigs to the reference genomes:"""
    metaSpades_assembly(args.name, read1, read2, args.assembly, args.isGT, args.num_processes)
    end_time = time.time()
    print("Time taken to assemble the contigs: ", end_time - start_time)

    """# **Get groundtruth repeats**"""
    start_time = time.time()
    if args.isGT:
        intra_Repeats, inter_Repeats = get_repeats(coords_file_name, myIDY=args.idy, myL=args.nL, myCN=args.cN)
    else:
        intra_Repeats = []
        inter_Repeats = []
    
    all_Repeats = intra_Repeats + inter_Repeats
    end_time = time.time()
    print("Time taken to find repeats: ", end_time - start_time)


    """# **Assembly Graph**"""
    start_time = time.time()
    G_hybrid= get_hybG(graphContigs, before_rr, name=args.name, num_processes=args.num_processes)
    end_time = time.time()
    print("Time taken to form and save assembly graph: ", end_time - start_time)

    """# **Feature extraction """
    start_time = time.time()
    feature_dict_G, df_unitig_node = get_nodeFeatures(G_hybrid, before_rr, args.num_processes)
    end_time = time.time()
    print("Time taken to features: ", end_time - start_time)

    """# **Save data**"""
    start_time = time.time()
    save_data(feature_dict_G, all_Repeats, df_unitig_node, args.name)
    end_time = time.time()
    print("Time taken to save data: ", end_time - start_time)

    sys.stdout = sys.__stdout__
