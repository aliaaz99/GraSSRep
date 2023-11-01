import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import time

def get_data_1(alignment1_name, alignment2_name):
    df1 = pd.read_table(alignment1_name, header=None,
                    names=["read", "flag", "unitig", "pos", "mapQ", "CIGAR", "unitigNext", "posNext", "templateL",
                           "SEQ"], usecols=[0, 2, 9])
    df2 = pd.read_table(alignment2_name, header=None,
                        names=["read", "flag", "unitig", "pos", "mapQ", "CIGAR", "unitigNext", "posNext", "templateL",
                            "SEQ"], usecols=[0, 2, 9])

    N = df1['SEQ'].first_valid_index()  # number of unitigs+2
    df1 = df1.iloc[N:, :]
    df2 = df2.iloc[N:, :]
    df1[['read', 'f/r']] = df1['read'].str.split("/", expand=True)
    df2[['read', 'f/r']] = df2['read'].str.split("/", expand=True)
    df1 = df1.drop(columns=['f/r', 'SEQ'])
    df2 = df2.drop(columns=['f/r', 'SEQ'])
    print("Number of unitigs: ", N - 2)
    print("Number of forward reads mapped: ", len(df1.read.unique()))
    print("Number of reverse reads mapped: ", len(df2.read.unique()))

    return df1, df2

def get_data_2(alignment1_name, alignment2_name):
    
    df1 = pd.read_table(alignment1_name, header=None, names=["read", "flag", "unitig"], usecols=[0, 2])
    df2 = pd.read_table(alignment2_name, header=None, names=["read", "flag", "unitig"], usecols=[0, 2])


    N = df1[~df1['read'].str.startswith('@')].index[0]

    df1 = df1.iloc[N:, :]
    df2 = df2.iloc[N:, :]
    df1[['f/r', 'read']] = df1['read'].str.split(".", expand=True)
    df2[['f/r', 'read']] = df2['read'].str.split(".", expand=True)
    df1 = df1.drop(columns=['f/r'])
    df2 = df2.drop(columns=['f/r'])
    print("Number of unitigs: ", N - 2)
    print("Number of forward reads mapped: ", len(df1.read.unique()))
    print("Number of reverse reads mapped: ", len(df2.read.unique()))

    return df1, df2

def get_repeats(coords_file_name, unitig_lengths, myIDY=100, myL=1, myCN=2):
    dfR = pd.read_csv(coords_file_name, skiprows=5, delimiter=' ', header=None, skipinitialspace=True)
    dfR.columns = ["S1", "E1", "r1", "S2", "E2", "r2", "Len1", "Len2", "r3", "IDY", "r4", "TAGS"]
    dfR = dfR.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    dfR[['TagR', 'TagU']] = dfR['TAGS'].str.split('\t', expand=True)
    dfR = dfR.drop(columns=['TAGS', 'r1', 'r2', 'r3', 'r4'])

    dfR.Len2 = dfR.Len2.astype(float)
    dfR.IDY = dfR.IDY.astype(float)
    dfR.TagU = dfR.TagU.astype(int)

    dfR['Len2_norm'] = dfR.apply(lambda row: row['Len2'] / unitig_lengths[row['TagU']]
                            if row['TagU'] in unitig_lengths else 0, axis=1)
          
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

def get_length(alignment_name):
    unitig_lengths = {}
    with open(alignment_name) as f:
        for line in f:
            if line.startswith('@SQ'):
                parts = line.strip().split('\t')
                unitig_id = int(parts[1][3:])
                length = int(parts[2][3:])
                unitig_lengths[unitig_id] = length

    return unitig_lengths

def get_coverage(alignment1_name, alignment2_name, unitig_lengths):

    coverage = {unitig_id: np.zeros((1, length)) for unitig_id, length in unitig_lengths.items()}

    # Process the forward reads
    with open(alignment1_name) as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            unitig_id = int(fields[2])
            pos = int(fields[3])
            read_len = len(fields[9])
            coverage[unitig_id][0][pos - 1:pos + read_len - 1] += 1
    # close the file
    f.close()
    del f

    # Process the reverse reads
    with open(alignment2_name) as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            unitig_id = int(fields[2])
            pos = int(fields[3])
            read_len = len(fields[9])
            coverage[unitig_id][0][pos - 1:pos + read_len - 1] += 1
    # close the file
    f.close()
    del f

    # Med and mean coverage
    mean_coverage = {unitig_id: np.mean(coverage[unitig_id]) for unitig_id in unitig_lengths.keys()}

    # get the maximum value in the dictionary and normalize
    max_coverage = max(mean_coverage.values())
    mean_coverage = {k: v / max_coverage for k, v in mean_coverage.items()}


    return mean_coverage

def get_RUmapping(df1, df2):
    group_dict1 = df1.groupby('read')['unitig'].unique().apply(list).to_dict()
    group_dict2 = df2.groupby('read')['unitig'].unique().apply(list).to_dict()

    return group_dict1, group_dict2

def get_repeatG(group_dict1, group_dict2):
    G_repeat = nx.Graph()

    group_dict1_repeat = {k: v for k, v in group_dict1.items() if len(v) > 1}
    group_dict2_repeat = {k: v for k, v in group_dict2.items() if len(v) > 1}

    for k, v in group_dict1_repeat.items():
        for i in range(len(v)):
            for j in range(i, len(v)):
                node1 = int(v[i])
                node2 = int(v[j])
                if node1 != node2:
                    edge = (node1, node2)
                    if G_repeat.has_edge(*edge):
                        # Edge already exists, update its weight
                        data = G_repeat.get_edge_data(*edge)
                        weight = data['weight'] + 1
                        G_repeat.add_edge(*edge, weight=weight)
                    else:
                        # Edge does not exist, add it with weight 1
                        G_repeat.add_edge(*edge, weight=1)

    for k, v in group_dict2_repeat.items():
        for i in range(len(v)):
            for j in range(i, len(v)):
                node1 = int(v[i])
                node2 = int(v[j])
                if node1 != node2:
                    edge = (node1, node2)
                    if G_repeat.has_edge(*edge):
                        # Edge already exists, update its weight
                        data = G_repeat.get_edge_data(*edge)
                        weight = data['weight'] + 1
                        G_repeat.add_edge(*edge, weight=weight)
                    else:
                        # Edge does not exist, add it with weight 1
                        G_repeat.add_edge(*edge, weight=1)

    print("G_repeat formed")

    return G_repeat

def get_adjG(group_dict1, group_dict2):
    G_adj = nx.Graph()
    # key1 is forward/reverse read name
    for key1 in group_dict1.keys():
        if key1 in group_dict2.keys():
            value1 = group_dict1[key1]  # name of unitigs forward read is mapped
            value2 = group_dict2[key1]  # name of unitigs reverse read is mapped
            for v1 in value1:
                for v2 in value2:
                    node1 = int(v1)
                    node2 = int(v2)
                    if node1 != node2:
                        edge = (node1, node2)
                        if G_adj.has_edge(*edge):
                            # Edge already exists, update its weight
                            data = G_adj.get_edge_data(*edge)
                            weight = data['weight'] + 1
                            G_adj.add_edge(*edge, weight=weight)
                        else:
                            # Edge does not exist, add it with weight 1
                            G_adj.add_edge(*edge, weight=1)

    print("G_adj formed")

    return G_adj

def get_hybG(G_repeat, G_adj, Repeat_U_0, Repeat_U_99, isPlot=1, isSave=1, name = 'shakya_1'):
    G_hybrid = nx.Graph()

    adj_edge_data = {edge: G_adj.get_edge_data(*edge) for edge in G_adj.edges}
    repeat_edge_data = {edge: G_repeat.get_edge_data(*edge) for edge in G_repeat.edges}

    common_edges = set(G_repeat.edges) & set(G_adj.edges)
    repeat_unique_edges = set(G_repeat.edges) - common_edges
    adj_unique_edges = set(G_adj.edges) - common_edges

    for edge in common_edges:
        weightT = adj_edge_data[edge]['weight'] + repeat_edge_data[edge]['weight']
        G_hybrid.add_edge(*edge, label='both', weight=weightT)

    for edge in repeat_unique_edges:
        G_hybrid.add_edge(*edge, label='repeat', weight=repeat_edge_data[edge]['weight'])

    for edge in adj_unique_edges:
        G_hybrid.add_edge(*edge, label='adj', weight=adj_edge_data[edge]['weight'])

    # print statictics of the weiths of the edges, mean, std, max, min
    weights = [d['weight'] for u, v, d in G_hybrid.edges(data=True)]
    print("Min weight of edges: ", np.min(weights))
    print("lower quartile weight of edges: ", np.percentile(weights, 25))
    print("Median weight of edges: ", np.median(weights))
    print("upper quartile weight of edges: ", np.percentile(weights, 75))
    print("Max weight of edges: ", np.max(weights))
    print("Mean weight of edges: ", np.mean(weights))
    

    # remove the edges with weight less than first quartile
    print("number of edges before removing edges: ", len(G_hybrid.edges))
    thresh_weight = np.percentile(weights, 25)
    G_hybrid.remove_edges_from([(u, v) for u, v, d in G_hybrid.edges(data=True) if d['weight'] < thresh_weight])
    print("number of edges after removing edges: ", len(G_hybrid.edges))

    if isPlot:
        # edge_colors = [("red" if d['label'] == 'repeat' else ("blue" if d['label'] == 'adj' else "black")) for u, v, d in
        #             G_hybrid.edges(data=True)]
        pos = nx.spring_layout(G_hybrid)
        node_colors = ['red' if n in Repeat_U_99 else ('yellow' if n in Repeat_U_0 else 'c') for n in G_hybrid.nodes()]
        plt.figure(figsize=(8, 8))
        nx.draw(G_hybrid, pos=pos, node_color=node_colors, node_size=300, width=1, font_size=8)
        title = "hybrid_graph"
        plt.title(title)
        red_patch = mpatches.Patch(color='red', label='inter-genome repeats')
        blue_patch = mpatches.Patch(color='yellow', label='intra-genome repeats')
        c_patch = mpatches.Patch(color='c', label='non-repeat')
        plt.legend(handles=[red_patch, blue_patch, c_patch], loc='upper right', fontsize=17.5)
        plt.savefig('Plots/' + name + '/' + title + '.png', dpi=300)

    if isSave:
        nx.write_edgelist(G_hybrid, 'Data/' + name + '/adj_matrix.coo', data=False, delimiter=' ')
        node_label_to_index = {label: index for index, label in enumerate(G_hybrid.nodes())}
        df_graph_map = pd.DataFrame(list(node_label_to_index.items()), columns=['unitig', 'idx'])
        df_graph_map.to_csv('Data/' + name + '/G_map.csv', index=False)


    print("G_hybrid formed and saved")

    return G_hybrid

def get_nodeFeatures(G_hybrid, alignment1_name, alignment2_name):
    node_label_to_index = {label: index for index, label in enumerate(G_hybrid.nodes())}
    df_unitig_node = pd.DataFrame(list(node_label_to_index.items()), columns=['unitig', 'idx'])

    """Simple degree":"""
    start_time = time.time()
    max_degree = max([G_hybrid.degree(n) for n in G_hybrid.nodes])
    print("max degree: ", max_degree)
    degree_dict_simple = {n: G_hybrid.degree(n) / max_degree for n in G_hybrid.nodes}
    print("Simple degree calculated")
    end_time = time.time()
    print("Time taken to calculate simple degree: ", end_time - start_time)

    """Weighted degree":"""
    start_time = time.time()
    total_edge_weight = sum([G_hybrid.get_edge_data(e[0], e[1])['weight'] for e in G_hybrid.edges])
    # get weighted degree of the nodes:
    degree_dict_weighted = {n: sum([G_hybrid.get_edge_data(n, m)['weight'] for m in G_hybrid.neighbors(n)]) / total_edge_weight for n in G_hybrid.nodes}
    print("Weighted degree calculated")
    end_time = time.time()
    print("Time taken to calculate weighted degree: ", end_time - start_time)

    """Unitig length:"""
    start_time = time.time()
    unitig_lengths = get_length(alignment1_name)
    print("Unitig lengths calculated")
    end_time = time.time()
    print("Time taken to calculate unitig lengths: ", end_time - start_time)

    """Coverage:"""
    start_time = time.time()
    mean_coverage = get_coverage(alignment1_name, alignment2_name, unitig_lengths)
    print("Coverage features calculated")
    end_time = time.time()
    print("Time taken to calculate coverage features: ", end_time - start_time)

    """Betweenness centrality:"""
    start_time = time.time()
    betweenness_dict = nx.betweenness_centrality(G_hybrid, weight=None)
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
    clustering_dict = dict(clustering_coefficient)
    print("Clustering coeff calculated")
    end_time = time.time()
    print("Time taken to calculate clustering coeff: ", end_time - start_time)

    print("All features calculated")
    feature_dict_G = {node:
                        [unitig_lengths[node],
                        mean_coverage[node],

                        degree_dict_simple[node],
                        degree_dict_weighted[node],
                        betweenness_dict[node],
                        core_dict[node],
                        clustering_dict[node],
                        ] for node in G_hybrid.nodes()}

    for unitig in unitig_lengths.keys():
        if unitig not in feature_dict_G.keys():
            new_row = {'unitig': unitig, 'idx': np.nan}
            df_unitig_node = df_unitig_node.append(new_row, ignore_index=True)
    
    return feature_dict_G, df_unitig_node

def save_data(feature_dict_G, intra_Repeats, inter_Repeats, all_Repeats, df_unitig_node, name = 'shakya_1'):
    X_G = np.array([feature_dict_G[k] for k in feature_dict_G])
    y_binary_G = np.array([1 if k in all_Repeats else 0 for k in feature_dict_G])
    y_3class_G = np.array([1 if k in intra_Repeats else (2 if k in inter_Repeats else 0) for k in feature_dict_G])
    N = X_G.shape[0]
    m = X_G.shape[1]

    """ save X and y for this graph:"""
    np.savez('Data/' + name + '/X_G.npz', X_G)
    np.savez('Data/' + name + '/y_binary.npz', y_binary_G)
    np.savez('Data/' + name + '/y_3class.npz', y_3class_G)

    print("X, y saved")
    print("Number of unitigs on graph: ", N)
    print("Number of features: ", m)
    print("Number of intra-repeat unitigs on the graph: ", np.sum(y_3class_G == 1))
    print("Number of inter-repeat unitigs on the graph: ", np.sum(y_3class_G == 2))
    print("Number of non-repeat unitigs on the graph: ", N - np.sum(y_binary_G))

    """ save X and y for non-graph unitigs:"""
    df_unitig_node['y_true_B'] = df_unitig_node.apply(lambda row: 1 if row['unitig'] in all_Repeats else 0, axis=1)
    df_unitig_node['y_true_nB'] = df_unitig_node.apply(lambda row: 1 if row['unitig'] in intra_Repeats else (2 if row['unitig'] in inter_Repeats else 0), axis=1)
    df_unitig_node.to_csv('Data/' + name + '/G_node_to_index.csv', index=False)
