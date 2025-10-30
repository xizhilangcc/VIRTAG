import  numpy as np
import  pickle as pkl
import  scipy.sparse as sp
import  torch
import pandas as pd
import os
import argparse


seed = 123
np.random.seed(seed)
torch.random.manual_seed(seed)

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--database_path', type=str, default = "./database")
parser.add_argument('--class_counter_path', type=str, default = "./database/class_counter.pkl")
parser.add_argument('--outpath', type=str, default = "result")
args = parser.parse_args()

adj = pkl.load(open(f"{args.outpath}/Cyber_data/contig.graph",'rb'))
indices_np = np.array(adj.nonzero())
indices = torch.tensor(indices_np, dtype=torch.long)
values = torch.tensor(adj.data, dtype=torch.float)
size = torch.Size(adj.shape)
adj = torch.sparse.FloatTensor(indices, values, size)
labels = pkl.load(open(f"{args.outpath}/Cyber_data/contig.label",'rb'))
labels = torch.tensor(labels)
labels = labels
train_mask = pkl.load(open(f"{args.outpath}/Cyber_data/contig.train_mask",'rb'))
train_mask = torch.tensor(train_mask)
test_mask = [x ^ 1 for x in train_mask]
test_mask = torch.tensor(test_mask)
zero_count = (test_mask == 1).sum().item()
test_dict = pkl.load(open(f"{args.outpath}/Cyber_data/contig.test_dict",'rb'))
class_counter = pkl.load(open(args.class_counter_path,'rb'))
class_counter = {value: key for key, value in class_counter.items()}


def label_propagation_accuracy(adj_matrix, labels, train_mask, test_mask, num_iterations):
    for iteration in range(num_iterations):
        test_predict = {}
        update_train_mask = train_mask.clone()
        update_labels = labels.clone()
        test_indices = torch.where(test_mask == 1)[0]
        number = 0
        for node in test_indices:
            neighbors = adj_matrix[node].to_dense().nonzero().squeeze()
            weights = adj_matrix[node].to_dense()[neighbors]
            train_neighbors = neighbors[train_mask[neighbors] == 1]
            train_weights = weights[train_mask[neighbors] == 1]
            if len(train_neighbors) == 0:
                continue
            update_train_mask[node] = 1
            number += 1
            neighbor_labels = labels[train_neighbors]
            unique_labels = neighbor_labels.unique()
            weighted_sums = {}
            for label in unique_labels:
                label_mask = neighbor_labels == label
                weighted_sums[label.item()] = train_weights[label_mask].sum().item() / label_mask.sum().item()
            predicted_label = max(weighted_sums, key=weighted_sums.get)
            test_sequence = test_dict[int(node)]
            predict = class_counter[predicted_label]
            if iteration == 0:
                test_predict[test_sequence] = predict
            else:
                test_predict[test_sequence] = predict+'_like'
            update_labels[node] = predicted_label

        update_test_mask = update_train_mask ^ 1
        df = pd.DataFrame(list(test_predict.items()), columns=["contig_name", "prediction"])
        csv_path = "result/pred/prediction.csv"
        if not os.path.exists(csv_path):
            df.to_csv(csv_path, index=False, mode="a", header=True)
        else:
            df.to_csv(csv_path, index=False, mode="a", header=False)

        labels = update_labels
        train_mask = update_train_mask
        test_mask = update_test_mask

    return train_mask, test_mask, labels

num_iterations = 10
update_train_mask, update_test_mask, update_labels = label_propagation_accuracy(
    adj, labels, train_mask, test_mask, num_iterations
)
