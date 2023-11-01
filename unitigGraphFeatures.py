from tools_graph import *
import sys
import time
import argparse as ap
import os
import json
 # """# **Initialize**"""

parser = ap.ArgumentParser()
parser.add_argument('--name', help='Folder containing reads and reference genomes (in Data folder)', default="shakya_1")
parser.add_argument('--idy', type=int, help='Identity for repeat detection in %', default=95)
parser.add_argument('--nL', type=float, help='Normalized length for repeat detection', default=0.95)
parser.add_argument('--cN', type=int, help='Copy number for repeat detection', default=2)
parser.add_argument('--isPlot', type=int, help='Plot the graph', default=1)
parser.add_argument('--isSave', type=int, help='Save the data', default=1)
args = parser.parse_args()

# name of input files:
alignment1_name = 'Data/' + args.name + '/alignment1.sam'
alignment2_name = 'Data/' + args.name + '/alignment2.sam'
coords_file_name = 'Data/' + args.name + '/unitig_mapping.coords'

# Creat folders to save plots and results of data
try:
    os.mkdir('Plots/' + args.name)
except:
    print("Plots folder already exists")

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
    """# **Read data**"""
    start_time = time.time()
    if args.name == "shakya_2":
        df1, df2 = get_data_2(alignment1_name, alignment2_name)
    else:
        df1, df2 = get_data_1(alignment1_name, alignment2_name)
    unitig_lengths = get_length(alignment1_name)
    all_unitig_lengths = [unitig_lengths[u] for u in unitig_lengths.keys()]
    intra_Repeats, inter_Repeats = get_repeats(coords_file_name, unitig_lengths=unitig_lengths, myIDY=args.idy, myL=args.nL, myCN=args.cN)
    all_Repeats = intra_Repeats + inter_Repeats
    end_time = time.time()
    print("Time taken to read data: ", end_time - start_time)

    """Set of unitigs for each read"""
    start_time = time.time()
    group_dict1, group_dict2 = get_RUmapping(df1, df2)
    end_time = time.time()
    print("Time taken to create group_dict: ", end_time - start_time)

    """# **Unitig Graph**"""
    start_time = time.time()
    G_repeat = get_repeatG(group_dict1, group_dict2)
    G_adj = get_adjG(group_dict1, group_dict2)
    G_hybrid= get_hybG(G_repeat, G_adj, intra_Repeats, inter_Repeats, isPlot=True, isSave=True, name=args.name)
    end_time = time.time()
    print("Time taken to form and save hybrid graph: ", end_time - start_time)

    """# **Feature extraction """
    start_time = time.time()
    feature_dict_G, df_unitig_node = get_nodeFeatures(G_hybrid, alignment1_name, alignment2_name)
    end_time = time.time()
    print("Time taken to features: ", end_time - start_time)

    """# **Save data**"""
    start_time = time.time()
    save_data(feature_dict_G, intra_Repeats, inter_Repeats, all_Repeats, df_unitig_node, args.name)
    end_time = time.time()
    print("Time taken to save data: ", end_time - start_time)

    sys.stdout = sys.__stdout__