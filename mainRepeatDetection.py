from toolsDetection import *
import torch
print('torch: ', torch.__version__)
from torch_geometric.data import Data
import numpy as np
from sklearn.metrics import confusion_matrix
import argparse as ap
import pandas as pd



# Parameters:
parser = ap.ArgumentParser()
parser.add_argument('--name', help='Folder containing reads and reference genomes (if available) in Data folder', default="shakya_1")
parser.add_argument('--p', type=int, help='p threshold to generate pseudo-labels', default=35)
parser.add_argument('--N_iter', type=int, help='Number of iteration', default=10)
parser.add_argument('--isSemiSupervised', type=int, help='Flag to use semi-supervised learning instead of self-supervised learning', default=0)
parser.add_argument('--noGNN', type=int, help='Flag to exclude GNN step or not', default=0)
parser.add_argument('--noRF', type=int, help='Flag to exclude RF step or not', default=0)
args = parser.parse_args()


name = args.name
N_tree = 100
thresholds = args.p
hidden_channel = 16
embSize = 8
N_epoch_GNN = 2001
gnnLr = 0.001

"""# **Data loading:**"""

print("Reading data: ")

# load the features and labels:
X_raw = np.load('Data/' + args.name + '/X_G.npz')
X_raw = X_raw['arr_0']
print("Number of features:", X_raw.shape[1])
y_binary = np.load('Data/' + args.name + '/y_binary.npz')
y_binary = y_binary['arr_0']
print("Number of repeats", np.sum(y_binary)) 
# load the graph data
contig_info = pd.read_csv('Data/' + args.name + '/G_node_to_index.csv')
g_map = pd.read_csv('Data/' + args.name + '/G_map.csv')
with open('Data/' + args.name + '/adj_matrix.txt', 'r') as file:
    lines = file.readlines()
edges = []
for line in lines:
    r, c = map(int, line.strip().split())
    edges.append((r, c))
df_edges = pd.DataFrame(edges, columns=['source', 'dest'])
df_edges['source_idx'] = df_edges['source'].map(g_map.set_index('contig')['idx'])
df_edges['dest_idx'] = df_edges['dest'].map(g_map.set_index('contig')['idx'])


# drop sequencing features from X_raw for the main self-supervised learnining method
X_seq = X_raw[:, 0:2]
if args.isSemiSupervised==0:
  X_raw = np.delete(X_raw, [0, 1, 6], axis=1)

# node attribute:
num_features = X_raw.shape[1]
X_attr = torch.tensor(X_raw, dtype=torch.float)
# create PyTorch tensor in COO format
edge_indices = np.array([df_edges['source_idx'], df_edges['dest_idx']])
print("edge_indices shape = ", edge_indices.shape)
edge_index = torch.tensor(edge_indices, dtype=torch.long)
# node labels
y_binary = torch.tensor(y_binary, dtype=torch.long)
# form the graph data
data_main = Data(x=X_attr, edge_index=edge_index, y_binary = y_binary)


"""**Repeat Detection:**"""
# print argparse parameters:
print("name:", args.name)
print("Testing threshold: ", args.p)
print("N_epoch_GNN:", N_epoch_GNN)
print("N_tree:", N_tree)
print("isSemiSupervised:", args.isSemiSupervised)
print("hidden_channel:", hidden_channel)
print("embSize:", embSize)
print("gnnLr:", gnnLr)
print("")


Test_Acc_gnn = []
Test_Sen_gnn = []
Test_Pre_gnn = []
Test_tpr_gnn = []
Test_fpr_gnn = []
Test_f1_gnn = []
length_per_gnn = []
conf_matrixes_gnn = []

Test_Acc_rf = []
Test_Sen_rf = []
Test_Pre_rf = []
Test_tpr_rf = []
Test_fpr_rf = []
Test_f1_rf= []
length_per_rf = []
conf_matrixes_rf= []
Test_Fi_rf = []

Test_Acc_Final_all = []
Test_Sen_Final_all = []
Test_Pre_Final_all = []
Test_tpr_Final_all = []
Test_fpr_Final_all = []
Test_f1_Final_all = []
length_per_Final_all = []
conf_matrixes_Final_all = []

for iter in range(args.N_iter):
    print("iter = ", iter)
    t = thresholds
    # initial labels based on length and coverage
    len_thresh_init = t
    coverage_thresh_init = 100 - t
    y_initial, train_mask = get_yinitial(X_seq, y_binary, len_thresh_init, coverage_thresh_init, args.isSemiSupervised)
    data_main.y_initial = torch.tensor(y_initial, dtype=torch.long)
    data_main.train_mask = train_mask


    model_GNN = GCN(in_channels=num_features, hidden_channels= hidden_channel, embSize=embSize, num_classes=2)
    X, y_true, y_pred_gnn, y_initial, confusion_matrix_gnn_iter, emb_all = embLearn(data_main, model_GNN, N_epoch_GNN, gnnLr, args.noGNN) 
    # performance of y_pred_gnn based on y_true
    conf_gnn = confusion_matrix(y_true, y_pred_gnn)
    acc_gnn = np.trace(conf_gnn) / np.sum(conf_gnn)
    sen_gnn = conf_gnn[1,1] / np.sum(conf_gnn[1,:]) if np.sum(conf_gnn[1,:])!=0 else 0
    pre_gnn = conf_gnn[1,1] / np.sum(conf_gnn[:,1]) if np.sum(conf_gnn[:,1])!=0 else 0
    tpr_gnn = conf_gnn[1,1] / np.sum(conf_gnn[1,:]) if np.sum(conf_gnn[1,:])!=0 else 0
    fpr_gnn = conf_gnn[0,1] / np.sum(conf_gnn[0,:]) if np.sum(conf_gnn[0,:])!=0 else 0
    f1_gnn = 2 * (pre_gnn * sen_gnn) / (pre_gnn + sen_gnn) if (pre_gnn + sen_gnn)!=0 else 0
    Test_Acc_gnn.append(acc_gnn)
    Test_Sen_gnn.append(sen_gnn)
    Test_Pre_gnn.append(pre_gnn)
    Test_tpr_gnn.append(tpr_gnn)
    Test_fpr_gnn.append(fpr_gnn)
    Test_f1_gnn.append(f1_gnn)
    conf_matrixes_gnn.append(conf_gnn)

    # train a random forest 
    if args.noRF==1:
       # assign  the gnn results to rf
        y_pred_rf = y_pred_gnn
        acc_rf = acc_gnn
        sen_rf = sen_gnn
        pre_rf = pre_gnn
        conf_rf = conf_gnn
        feature_importance_rf = np.zeros(X.shape[1])
    else:
      acc_rf, sen_rf, pre_rf, conf_rf, y_pred_rf, feature_importance_rf = rfLearn(X, y_initial, y_true, train_mask, N_tree)

    tpr_rf = sen_rf
    fpr_rf = conf_rf[0,1] / np.sum(conf_rf[0,:]) if np.sum(conf_rf[0,:])!=0 else 0
    f1_rf = 2 * (pre_rf * sen_rf) / (pre_rf + sen_rf) if (pre_rf + sen_rf)!=0 else 0
    Test_Acc_rf.append(acc_rf)
    Test_Sen_rf.append(sen_rf)
    Test_Pre_rf.append(pre_rf)
    Test_tpr_rf.append(tpr_rf)
    Test_fpr_rf.append(fpr_rf)
    Test_f1_rf.append(f1_rf)
    conf_matrixes_rf.append(conf_rf)
    Test_Fi_rf.append(feature_importance_rf)

    if args.isSemiSupervised==1:
      y_final = y_pred_rf
    else:
      # get index of the points predicted 0 and 1
      indeces_0 = np.where(y_pred_rf == 0)[0]
      indeces_1 = np.where(y_pred_rf == 1)[0]

      # find the threshold of length and coverage of the points predicted 0
      thresh_len_0 = np.percentile(X_seq[indeces_0,0], len_thresh_init)
      thresh_mC_0 = np.percentile(X_seq[indeces_0,1], coverage_thresh_init)

      # find the threshold of length and coverage of the points predicted 1
      thresh_len_1 = np.percentile(X_seq[indeces_1,0], 100 - len_thresh_init)
      thresh_mC_1 = np.percentile(X_seq[indeces_1,1], coverage_thresh_init) 

      # update the labels based on the upper quartile of length and coverage of the points predicted 0 by GNN
      y_final = np.zeros(len(y_true))
      for i in range(len(y_true)):
          len_i = X_seq[i,0]
          mC_i = X_seq[i,1]
          if i in train_mask: # keep training labels
            y_final[i] = y_initial[i]
          else:
            if y_pred_rf[i]==0:
                if len_i < thresh_len_0 and mC_i > thresh_mC_0:
                  y_final[i] = 1
                else:
                  y_final[i] = 0
            else:
                if len_i > thresh_len_1 and mC_i < thresh_mC_1:
                    y_final[i] = 0
                else:
                    y_final[i] = 1


    contig_info['y_pred_B'] = contig_info.apply(lambda row: 1 if ((not pd.isna(row['idx'])) and (y_final[int(row['idx'])])) else 0, axis=1)      
    contig_info.loc[contig_info['contig'].isna(), 'y_pred_B'] = 0

    contig_info['train_mask'] = contig_info.apply(lambda row: 1 if row['idx'] in train_mask else 0, axis=1)

    # final evaluation on all of the contigs using y_true_B and y_pred_B column of contig_info
    y_true_B = np.array(contig_info['y_true_B'])
    y_pred_B = np.array(contig_info['y_pred_B'])

    conf_final_all_u = confusion_matrix(y_true_B, y_pred_B)
    acc_final_all_u = np.trace(conf_final_all_u) / np.sum(conf_final_all_u)
    sen_final_all_u = conf_final_all_u[1,1] / np.sum(conf_final_all_u[1,:]) if np.sum(conf_final_all_u[1,:])!=0 else 0
    tpr_final_all_u = conf_final_all_u[1,1] / np.sum(conf_final_all_u[1,:]) if np.sum(conf_final_all_u[1,:])!=0 else 0
    fpr_final_all_u = conf_final_all_u[0,1] / np.sum(conf_final_all_u[0,:]) if np.sum(conf_final_all_u[0,:])!=0 else 0
    pre_final_all_u = conf_final_all_u[1,1] / np.sum(conf_final_all_u[:,1]) if np.sum(conf_final_all_u[:,1])!=0 else 0
    f1_final_all_u = 2 * (pre_final_all_u * sen_final_all_u) / (pre_final_all_u + sen_final_all_u) if (pre_final_all_u + sen_final_all_u)!=0 else 0

    Test_Acc_Final_all.append(acc_final_all_u)
    Test_Sen_Final_all.append(sen_final_all_u)
    Test_Pre_Final_all.append(pre_final_all_u)
    Test_tpr_Final_all.append(tpr_final_all_u)
    Test_fpr_Final_all.append(fpr_final_all_u)
    Test_f1_Final_all.append(f1_final_all_u)
    conf_matrixes_Final_all.append(conf_final_all_u)
    

    if f1_final_all_u>= np.max(np.array(Test_f1_Final_all)):
      print("Best f1 till now:", f1_final_all_u)
      # update contig_info with adding predicted labels and trainmask as new column based on the 'idx' column of the dataframe
      contig_info.to_csv('Data/' + args.name + '/final_pred.csv', index=False)


print("")
print("Performance of initial labels:")
# performance of y_initial based on y_binary
conf_matrix_init = confusion_matrix(y_binary, y_initial)
acc_init = 100 * np.trace(conf_matrix_init) / np.sum(conf_matrix_init)
sen_init = 100 * conf_matrix_init[1,1] / np.sum(conf_matrix_init[1,:]) if np.sum(conf_matrix_init[1,:])!=0 else 0
pre_init = 100 * conf_matrix_init[1,1] / np.sum(conf_matrix_init[:,1]) if np.sum(conf_matrix_init[:,1])!=0 else 0
FPR_init = 100 * conf_matrix_init[0,1] / np.sum(conf_matrix_init[0,:]) if np.sum(conf_matrix_init[0,:])!=0 else 0
f1_init = 2 * (pre_init * sen_init) / (pre_init + sen_init) if (pre_init + sen_init)!=0 else 0
print("Initial labels (Baseline) performance on whole data:")
print("Accuracy:", acc_init)
print("Sensitivity:", sen_init)
print("Precision:", pre_init)
print("TPR:", sen_init)
print("FPR:", FPR_init)
print("F1:", f1_init)
print("Confusion Matrix: \n", conf_matrix_init)

print("")
print("Performance of GNN:")
print("Accuracy:", np.mean(np.array(Test_Acc_gnn))  * 100)
print("Sensitivity:",  np.mean(np.array(Test_Sen_gnn))  * 100)
print("Precision:",  np.mean(np.array(Test_Pre_gnn))  * 100)
print("TPR:",  np.mean(np.array(Test_tpr_gnn))  * 100)
print("FPR:",  np.mean(np.array(Test_fpr_gnn))  * 100)
print("F1:",  np.mean(np.array(Test_f1_gnn))  * 100)
print("Confusion Matrix: \n",  np.mean(np.array(conf_matrixes_gnn), axis=0))

print("")
print("Performance of RF:")
print("Accuracy:", np.mean(np.array(Test_Acc_rf))  * 100)
print("Sensitivity:",  np.mean(np.array(Test_Sen_rf))  * 100)
print("Precision:",  np.mean(np.array(Test_Pre_rf))  * 100)
print("TPR:",  np.mean(np.array(Test_tpr_rf))  * 100)
print("FPR:",  np.mean(np.array(Test_fpr_rf))  * 100)
print("F1:",  np.mean(np.array(Test_f1_rf))  * 100)
print("Confusion Matrix: \n",  np.mean(np.array(conf_matrixes_rf), axis=0))
print("Feature Importance: \n",  np.mean(np.array(Test_Fi_rf), axis=0))

print("")
print("All Final Results on all contigs after fine tuning:")
print("Accuracy:", np.mean(np.array(Test_Acc_Final_all))  * 100)
print("Sensitivity:",  np.mean(np.array(Test_Sen_Final_all))  * 100)
print("Precision:",  np.mean(np.array(Test_Pre_Final_all))  * 100)
print("TPR:",  np.mean(np.array(Test_tpr_Final_all))  * 100)
print("FPR:",  np.mean(np.array(Test_fpr_Final_all))  * 100)
print("F1:",  np.mean(np.array(Test_f1_Final_all))  * 100)
print("Confusion Matrix: \n",  np.mean(np.array(conf_matrixes_Final_all), axis=0))

