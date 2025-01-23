import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
import pandas as pd
import scipy


class dataset(Dataset):
    """\
        Description:
        ------------
            Create Pytorch Dataset

        Parameters:
        ------------
            counts: gene count. Type: numpy ndarrary
            anchor: anchor index. Type: numpy ndarray.
        
        Return:
        ------------
            Dataset
        """
    def __init__(self, counts, anchor = None):

        # assert not len(counts) == 0, "Count is empty"
        B = counts > 0
        indexes = scipy.sparse.find(B)
        self.indices = torch.tensor([indexes[0], indexes[1]], dtype=torch.long)
        self.values = torch.tensor(counts.data, dtype=torch.float32)
        self.shape = torch.Size(counts.shape)
        self.counts = torch.sparse.FloatTensor(self.indices, self.values, self.shape)
        # self.counts = torch.FloatTensor(counts)

        self.is_anchor = np.zeros(self.counts.shape[0]).astype("bool")    
        if anchor is not None:
            self.is_anchor[anchor] = True
        
        self.is_anchor = torch.tensor(self.is_anchor)
                   
    def __len__(self):
        return self.counts.shape[0]
    
    def __getitem__(self, idx):
        # data original data, index the index of cell, label, corresponding labels, batch, corresponding batch number
        row_index = idx  # 要提取的行索引

        row_indices = self.indices[0, :]  # 索引中的行索引
        col_indices = self.indices[1, :]
        row_idx = row_indices[row_indices == idx]
        col_idx = col_indices[row_indices == idx]
        index = self.indices[:,row_indices == idx]
        row_values = self.values[row_indices == row_index]
        dense_matrix = torch.sparse.FloatTensor(index, row_values, self.shape).to_dense()
        matrix = dense_matrix[idx,:]
        sample = {"count": matrix, "index": idx, "is_anchor": self.is_anchor[idx]}
        # sample = {"count": self.counts[idx,:], "index": idx, "is_anchor": self.is_anchor[idx]}
        return sample
