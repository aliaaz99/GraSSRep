import numpy as np
import pandas as pd
import pyfastg
import networkx as nx
from Bio import SeqIO
import time
import subprocess
import multiprocessing
import itertools
from multiprocessing import Pool


def metaSpades_assembly(folder_name, read1, read2, isAssembly, isGT, num_processes):
    genome_name = 'Data/' + folder_name + '/ref_genome.fasta'
    spades_out = 'Data/' + folder_name  + '/assemblies'
    myUnitig = spades_out + '/before_rr.fasta'
    nucmer_name = spades_out + '/contig_mapping'
    nucmer_repeat = nucmer_name + '.coords'
    k = 55
    
    if isAssembly:
        spades_cmd = ['spades.py' + ' -1 ' + read1 + ' -2 ' + read2 + ' -o ' + spades_out + ' -t ' + str(num_processes) + ' --meta' + ' --disable-rr ' + ' -m 70']
        with open('Data/' + folder_name + '/spades_terminal_file.txt', "w") as outfile:
            p_spades = subprocess.run(spades_cmd, shell=True, stdout=outfile, stderr=outfile)

    if isGT:
        command_nucmer_1 = ['nucmer' + ' --maxmatch' + ' --nosimplify' + ' -l ' + str(int(k)-1) + ' -c ' + str(int(k)-1) + ' --prefix=' + nucmer_name + ' ' + genome_name + ' ' + myUnitig]
        command_nucmer_2 = ['show-coords' + ' -r ' + nucmer_name + '.delta' + ' > ' + nucmer_repeat]
        with open('Data/' + folder_name + '/nucmer_terminal_file1.txt', "w") as outfile:
            p_nucmer1 = subprocess.run(command_nucmer_1, shell=True, stdout=outfile, stderr=outfile)

        with open('Data/' + folder_name + '/nucmer_terminal_file2.txt', "w") as outfile:
            p_nucmer2 = subprocess.run(command_nucmer_2, shell=True, stdout=outfile, stderr=outfile)  

    print("Assembly done")



def get_repeats(coords_file_name, myIDY=100, myL=100, myCN=2):
    dfR = pd.read_csv(coords_file_name, skiprows=5, delimiter=' ', header=None, skipinitialspace=True)
    dfR.columns = ["S1", "E1", "r1", "S2", "E2", "r2", "Len1", "Len2", "r3", "IDY", "r4", "TAGS"]
    dfR = dfR.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    dfR[['TagR', 'TagU']] = dfR['TAGS'].str.split('\t', expand=True)
    dfR[['NODE', 'TagU', 'length', 'LenU', 'coverage','coverageU']] = dfR['TagU'].str.split('_', expand=True)
    dfR = dfR.drop(columns=['S1', 'E1', 'r1', 'S2', 'E2', 'r2', 'Len1', 'TAGS', 'r3', 'r4', 'NODE', 'length', 'coverage', 'coverageU'])

    dfR.Len2 = dfR.Len2.astype(int)
    dfR.IDY = dfR.IDY.astype(float)
    dfR.TagU = dfR.TagU.astype(int)
    dfR.LenU = dfR.LenU.astype(int)
    dfR['Len2_norm'] = dfR['Len2'] / dfR['LenU'] * 100
    
    filtered_rows = dfR[dfR['TagU'].duplicated(keep=False) & (dfR['IDY'] >= myIDY) & (dfR['Len2_norm'] >= myL)]
    tagu_counts = filtered_rows.groupby('TagU').size()
    filtered_df = filtered_rows[filtered_rows['TagU'].isin(tagu_counts[tagu_counts >= myCN].index)]

    Repeat_U_0 = []    
    Repeat_U_99 = []  
    for tagu, group in filtered_df.groupby('TagU'):
        tagr_values = group['TagR'].unique()
        if len(tagr_values) > 1:
            Repeat_U_99.append(tagu)
        else:
            Repeat_U_0.append(tagu)


    print("Number of Intra Repeat unitigs:", len(Repeat_U_0))
    print("Number of Inter Repeat unitigs:", len(Repeat_U_99))
    return Repeat_U_0, Repeat_U_99

def find_map_chunck(before_rr_dict):
    node_id = []
    contig_id = []
    for contig_name, sequence1 in before_rr_dict.items():
        for node_name, sequence2 in fastg_dict.items():
            if sequence2 == sequence1:
                node_name_map = node_name
                break
        
        name_1_split = contig_name.split('_')
        name_2_split = node_name_map.split('_')

        contig_id_i = name_1_split[1]
        node_id_i = name_2_split[1]
        
        node_id.append(node_id_i+'+')
        contig_id.append(int(contig_id_i))
        node_id.append(node_id_i+'-')
        contig_id.append(int(contig_id_i))
    
    return node_id, contig_id

def map_sequence_names(before_rr, fastg, name, num_processes):
    global fastg_dict
    before_rr_dict = {}
    fastg_dict = {}
    for record in SeqIO.parse(before_rr, "fasta"):
        before_rr_dict[record.id] = record.seq
    for record in SeqIO.parse(fastg, "fasta"):
        fastg_dict[record.id] = record.seq

    print("Len of before_rr_dict: ", len(before_rr_dict))
    print("Len of fastg_dict: ", len(fastg_dict))

    # split the before_rr dictionary to num_processes smaller dictionaries:
    print("num_processes: ", num_processes)
    smaller_dicts = [{} for _ in range(num_processes)]
    print("smaller_dicts len: ", len(smaller_dicts))
    k = len(before_rr_dict) // num_processes + 1
    print("k: ", k)
    for i, (key, value) in enumerate(before_rr_dict.items()):
        smaller_dicts[i // k][key] = value
            
    args_list = []
    for i in range(num_processes):
        args_list.append(smaller_dicts[i])
    print("arg_list len: ", len(args_list))

    node_id = []
    contig_id = []
    with multiprocessing.Pool(processes=num_processes) as p:
        results = p.map(find_map_chunck, args_list)
        for node_id_i, contig_id_i in results:
            node_id.extend(node_id_i)
            contig_id.extend(contig_id_i)

    # make a dictionary from node_id to contig_id:
    mapping_dict = dict(zip(node_id, contig_id))
    mapping_dict_df = pd.DataFrame(list(mapping_dict.items()), columns=['bandage', 'contig'])
    mapping_dict_df.to_csv('Data/' + name + '/Bandage_map.csv', index=False)
    return mapping_dict

def get_hybG(graphContigs, before_rr, name = 'shakya_1', num_processes=multiprocessing.cpu_count()):
    G_hybrid = pyfastg.parse_fastg(graphContigs)
    mapping_dict = map_sequence_names(before_rr, graphContigs, name, num_processes)
    G_hybrid = nx.relabel_nodes(G_hybrid, mapping_dict)
    G_hybrid = G_hybrid.to_undirected()
    nx.write_edgelist(G_hybrid, 'Data/' + name + '/adj_matrix.txt', data=False, delimiter=' ')
    node_label_to_index = {label: index for index, label in enumerate(G_hybrid.nodes())}
    df_graph_map = pd.DataFrame(list(node_label_to_index.items()), columns=['contig', 'idx'])
    df_graph_map.to_csv('Data/' + name + '/G_map.csv', index=False)
    print("G_hybrid formed and saved")
    print(G_hybrid)
    return G_hybrid

def get_length_cov(before_rr):
    unitig_lengths = {}
    unitig_coverages = {}
    with open(before_rr) as f:
        for line in f:
            if line.startswith('>NODE'):
                parts = line.strip().split('_')
                unitig_id = int(parts[1])
                length = int(parts[3])
                coverage = float(parts[5])
                unitig_lengths[unitig_id] = length
                unitig_coverages[unitig_id] = coverage

    return unitig_lengths, unitig_coverages

def chunks(l, n):
    """Divide a list of nodes `l` in `n` chunks"""
    l_c = iter(l)
    while 1:
        x = tuple(itertools.islice(l_c, n))
        if not x:
            return
        yield x

def betweenness_centrality_parallel(G, processes=None):
    """Parallel betweenness centrality  function"""
    p = Pool(processes=processes)
    node_divisor = len(p._pool) * 4
    node_chunks = list(chunks(G.nodes(), G.order() // node_divisor))
    num_chunks = len(node_chunks)
    bt_sc = p.starmap(
        nx.betweenness_centrality_subset,
        zip(
            [G] * num_chunks,
            node_chunks,
            [list(G)] * num_chunks,
            [True] * num_chunks,
            [None] * num_chunks,
        ),
    )

    # Reduce the partial solutions
    bt_c = bt_sc[0]
    for bt in bt_sc[1:]:
        for n in bt:
            bt_c[n] += bt[n]
    return bt_c

def get_nodeFeatures(G_hybrid, before_rr, num_processes=multiprocessing.cpu_count()):
    
    """Unitig length and coverage:"""
    start_time = time.time()
    unitig_lengths, unitig_coverages = get_length_cov(before_rr)  
    end_time = time.time()  
    print("Time taken to calculate length and coverage features: ", end_time - start_time)

    """Simple degree":"""
    start_time = time.time()
    max_degree = max([G_hybrid.degree(n) for n in G_hybrid.nodes])
    print("max degree: ", max_degree)
    degree_dict_simple = {n: G_hybrid.degree(n) / max_degree for n in G_hybrid.nodes}
    print("Simple degree calculated")
    end_time = time.time()
    print("Time taken to calculate simple degree: ", end_time - start_time)

    # remove the self loops:
    G_hybrid.remove_edges_from(nx.selfloop_edges(G_hybrid))
    
    """Betweenness centrality:"""
    start_time = time.time()
    betweenness_dict = betweenness_centrality_parallel(G_hybrid, num_processes)
    betweenness_dict = {node: betweenness_dict[node] for node in G_hybrid.nodes()}
    print("Betweenness centrality calculated")
    end_time = time.time()
    print("Time taken to calculate betweenness centrality: ", end_time - start_time)

    """K-Core number:"""
    start_time = time.time()
    core_dict = nx.core_number(G_hybrid)
    core_dict = {node: core_dict[node] for node in G_hybrid.nodes()}
    print("K-Core number calculated")
    end_time = time.time()
    print("Time taken to calculate K-Core number: ", end_time - start_time)

    """Clustering coeff:"""
    start_time = time.time()
    clustering_coefficient = nx.clustering(G_hybrid, weight=None)
    clustering_dict = {node: clustering_coefficient[node] for node in G_hybrid.nodes()}
    print("Clustering coeff calculated")
    end_time = time.time()
    print("Time taken to calculate clustering coeff: ", end_time - start_time)

    """skewed links:"""
    start_time = time.time()
    skewed_edges = {}
    for node in G_hybrid.nodes():
        num_neighs = len(list(G_hybrid.neighbors(node)))
        if num_neighs == 0:
            skewed_edges[node] = 0
        else:
            s_count = 0
            for neighs in G_hybrid.neighbors(node):
                if unitig_coverages[node] >= 2*unitig_coverages[neighs]:
                        s_count += 1
            skewed_edges[node] = s_count/ num_neighs
    print("Skewed links calculated")
    end_time = time.time()
    print("Time taken to calculate skewed links: ", end_time - start_time)

    print("All features calculated")
    feature_dict_G = {node:
                        [unitig_lengths[node],
                        unitig_coverages[node],

                        degree_dict_simple[node],
                        betweenness_dict[node],
                        core_dict[node],
                        clustering_dict[node],
                        skewed_edges[node]
                        ] for node in G_hybrid.nodes()}

    node_label_to_index = {node: index for index, node in enumerate(G_hybrid.nodes())}
    df_unitig_node = pd.DataFrame(list(node_label_to_index.items()), columns=['contig', 'idx'])
    
    return feature_dict_G, df_unitig_node

def save_data(feature_dict_G, all_Repeats, df_unitig_node, name = 'shakya_1'):
    X_G = np.array([feature_dict_G[k] for k in feature_dict_G])
    y_binary_G = np.array([1 if k in all_Repeats else 0 for k in feature_dict_G])
    N = X_G.shape[0]
    m = X_G.shape[1]

    """ save X and y for this graph:"""
    np.savez('Data/' + name + '/X_G.npz', X_G)
    np.savez('Data/' + name + '/y_binary.npz', y_binary_G)
    print("X, y saved")

    print("Number of unitigs on graph: ", N)
    print("Number of features: ", m)
    print("Number of non-repeat unitigs on the graph: ", N - np.sum(y_binary_G))

    """ save X and y for non-graph unitigs:"""
    df_unitig_node['y_true_B'] = df_unitig_node.apply(lambda row: 1 if row['contig'] in all_Repeats else 0, axis=1)
    df_unitig_node.to_csv('Data/' + name + '/G_node_to_index.csv', index=False)
