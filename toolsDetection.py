import numpy as np
import torch
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv
from sklearn.metrics import confusion_matrix, accuracy_score, recall_score, precision_score
from sklearn.ensemble import RandomForestClassifier
from torch_geometric.loader import NodeLoader

    
class GCN(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, embSize, num_classes):
        super(GCN, self).__init__()
        self.conv1 = SAGEConv(in_channels, hidden_channels)
        self.bn1 = torch.nn.BatchNorm1d(hidden_channels)
        self.conv2 = SAGEConv(hidden_channels, embSize)
        self.bn2 = torch.nn.BatchNorm1d(embSize)
        self.fc = torch.nn.Linear(embSize, num_classes)
        
    def forward(self, x, edge_index):
        x0 = self.conv1(x, edge_index)
        x0 = self.bn1(x0)
        x0 = F.relu(x0)

        x0 = self.conv2(x0, edge_index)
        x0 = self.bn2(x0)
        x0 = F.relu(x0)

        x1 = self.fc(x0)

        x1 = F.softmax(x1, dim=1)

        return x0, x1
    
def get_yinitial(X, y_real, len_thresh, coverage_thresh, isRandom=0):
    flags_ones = np.zeros(X.shape[0])
    flags_zeros = np.zeros(X.shape[0])
    # length flag:
    lower_quartile_len = np.percentile(X[:, 0], len_thresh)
    ones_indices_len = X[:, 0] < lower_quartile_len
    flags_ones[ones_indices_len] += 1

    upper_quartile_len = np.percentile(X[:, 0], 100 - len_thresh)
    zeros_indices_len = X[:, 0] > upper_quartile_len
    flags_zeros[zeros_indices_len] += 1

    # mean coverage flag:
    upper_quartile_meanC = np.percentile(X[:, 1], coverage_thresh)
    ones_indices_mC = X[:, 1] > upper_quartile_meanC
    flags_ones[ones_indices_mC] += 1

    lower_quartile_meanC = np.percentile(X[:, 1], coverage_thresh)
    zeros_indices_mC = X[:, 1] < lower_quartile_meanC
    flags_zeros[zeros_indices_mC] += 1

    # Find the datapoints that are flagged by both of the features
    labelled_indices_ones = np.where(flags_ones >= 2)[0]
    labelled_indices_zeros = np.where(flags_zeros >= 2)[0]
    # training indeces:
    train_indeces = labelled_indices_ones.tolist() + labelled_indices_zeros.tolist()    

    if isRandom:
        N_train = int(len(train_indeces))
        train_indeces = np.random.choice(X.shape[0], N_train, replace=False)
        train_mask = np.zeros(X.shape[0], dtype=bool)
        train_mask[train_indeces] = True
        train_mask = torch.tensor(train_mask, dtype=torch.bool)
        y_initial = y_real
    else:
        train_mask = np.zeros(X.shape[0], dtype=bool)
        train_mask[train_indeces] = True
        train_mask = torch.tensor(train_mask, dtype=torch.bool)

        # initial labels:
        y_initial = np.zeros(X.shape[0], dtype=bool)
        y_initial[labelled_indices_ones] = True
        y_initial = y_initial.astype(int)
        # precision of ones:
        ones_mask = np.zeros(X.shape[0], dtype=bool)
        ones_mask[labelled_indices_ones] = True
        print("traning data precision of ones:", 100 * precision_score(y_real, ones_mask))
        # precision of zeros
        zeros_mask = np.ones(X.shape[0], dtype=bool)
        zeros_mask[labelled_indices_zeros] = False
        conf_zeros = confusion_matrix(y_real, zeros_mask)
        zeros_precision = conf_zeros[0,0] / (conf_zeros[0,0] + conf_zeros[1,0])
        print("traning data precision of zeros:", 100 *  zeros_precision)


    return y_initial, train_mask

def embLearn(data_in, model_GNN, N_epoch_GNN, gnnLr, noGNN):
    print("learning embeddings...")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Using device:", device)

    x_all, edge_index_all, y_all, y_initial, train_mask = (
    data_in.x.to(device),
    data_in.edge_index.to(device),
    data_in.y_binary.to(device),
    data_in.y_initial.to(device),
    data_in.train_mask.to(device),)
    
    model_GNN = model_GNN.to(device)

    optimizer_GNN = torch.optim.Adam(model_GNN.parameters(), lr=gnnLr)

    criterion_GNN = torch.nn.CrossEntropyLoss()

    N_epoch_GNN_set = 1 if noGNN else N_epoch_GNN
    model_GNN.train()
    for epoch in range(N_epoch_GNN_set):
        optimizer_GNN.zero_grad()
        emb_, out_GNN = model_GNN(x_all, edge_index_all)
        loss = criterion_GNN(out_GNN[train_mask], y_initial[train_mask])
        loss.backward()
        optimizer_GNN.step()
        epoch_pred = out_GNN[train_mask].argmax(dim=1).cpu().numpy()  # Move predictions back to CPU
        epoch_real = y_initial[train_mask].cpu().numpy()  # Move labels back to CPU
        loss_epoch = loss.item()
        
        if epoch % 100 == 0:
            acc = accuracy_score(epoch_real, epoch_pred)
            rec = recall_score(epoch_real, epoch_pred, average='macro')
            print('y_initial: Epoch: {:03d}, loss: {:.5f}, Acc: {:.5f}, Rec: {:.5f}'.format(epoch, loss_epoch, acc, rec))

    model_GNN.eval()
    with torch.no_grad():
        emb_GNN_all, out = model_GNN(x_all, edge_index_all)
        y_pred = out.argmax(dim=1).cpu().numpy()  # Move predictions back to CPU
        y_true = y_all.cpu().numpy()  # Move labels back to CPU

    confusion_matrix_gnn = confusion_matrix(y_true, y_pred)  # Convert predictions to CPU for the confusion matrix
    emb_GNN_all = emb_GNN_all.cpu().numpy()  # Move embeddings back to CPU

    # Create X_train and X_test by combining node features and GNN embeddings
    x_features = x_all.cpu().numpy()
    y_initial = y_initial.cpu().numpy()

    if noGNN:
        X = x_features
    else:
        X = np.concatenate((x_features, emb_GNN_all), axis=1)


    return X, y_true, y_pred, y_initial, confusion_matrix_gnn, emb_GNN_all

def rfLearn(X, y_initial, y_true, train_mask, N_tree):
    X_train = X[train_mask]
    y_train = y_initial[train_mask]
    rf = RandomForestClassifier(n_estimators=N_tree, n_jobs = 30)
    rf.fit(X_train, y_train)

    y_pred_rf = rf.predict(X)
    conf_rf = confusion_matrix(y_true, y_pred_rf)
    acc_rf = np.trace(conf_rf) / np.sum(conf_rf) 
    sen_rf = conf_rf[1,1] / np.sum(conf_rf[1,:]) if np.sum(conf_rf[1,:])!=0 else 0
    pre_rf = conf_rf[1,1] / np.sum(conf_rf[:,1]) if np.sum(conf_rf[:,1])!=0 else 0
    
    return acc_rf, sen_rf, pre_rf, conf_rf, y_pred_rf, rf.feature_importances_
