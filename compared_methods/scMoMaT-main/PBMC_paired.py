#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys, os

import numpy as np
from umap import UMAP
import time
import torch
import matplotlib.pyplot as plt
import pandas as pd  
import scipy.sparse as sp
import scanpy as sc

import scmomat 

plt.rcParams["font.size"] = 10


# ## Step 1: Load dataset
# Make sure that the dataset get stored in a dictionary (`counts` in the script), with the item:
# * `nbatches`: total number of batches
# * `feats_name`: a dictionary storing the feature names for different modalities, e.g. `{"rna": np.array([aaaa, bbbb,...]), "atac": np.array([ccc, ddd, ...]), "protein": np.array([eee, fff, ...])}`
# * `modality 1` (e.g. `rna` in the example below): a `list` store all data matrices correspond to modality 1, a element correspond to one batch, elements are ordered in the list following the ordering of batches. **The batches that does not have count matrix in corresponding modality is inserted `None` as placeholder**
# * `modality 2` (e.g. `atac` in the example below): requirement the same as above.
# * `modality 3` (e.g. `protein`): requirement the same as above.
# * ...
# 
# #### Note:
# * The number of item in the `feats_name` should match the number of modalities in `counts`. 
# * The length of lists in `modality 1`, `modality 2`, `modality 3`, etc should have the same length, which is equal to `nbatches`. (missing matrices are denoted as `None` as explained above). **The matrices must be ordered according to the batch ID in each list,** so that scMoMaT can detact parallel sequenced batches.
# * The data matrix in each modality (each list) should have the same set of features. You can do thie by 1. using the intersection of all genes/proteins in different batches; 2. Remap the chromatin regions according to the peaks of one batch. 
# * The data matrix is of the shape `(ncells, nfeatures)`, and preprocessed with `utils.preprocess()` function.
# 
# One example template is shown as below, note that the features of data matrices are matched in advance. The connection matrix (`GxR.npz`) that is used to generate pseudo-count is also provided in advance (used the code in folder `calc_pseudo_count/calc_pseudo_count.R`). Feel free to modify on the template to use your dataset.
# 

# In[2]:


# data_dir = "./data/real/MOp_5batches/"

# # number of batches
# n_batches = 5

# # obtain the feature name
# genes = pd.read_csv(data_dir + "genes.txt", header = None).values.squeeze()
# regions = pd.read_csv(data_dir + "regions.txt", header = None).values.squeeze()
# feats_name = {"rna": genes, "atac": regions}

data_id = 'PBMC_paired'
rna = sc.read_h5ad('../../data/'+ data_id +'/GASM/rna.h5ad')
atac = sc.read_h5ad('../../data/'+ data_id +'/GASM/atac.h5ad')
n_batches = 2

reg = sc.read_h5ad('../datasets/scDART/data/' + data_id + '.h5ad')
A = np.array(reg.X.todense())

# obtain the feature name
genes = np.array(reg.obs_names)
regions = np.array(reg.var_names)

feats_name = {"rna": genes, "atac": regions}

# READ IN THE COUNT MATRICES
# subsample the original dataset to save the training time
subsample = 10
# scATAC-seq of batch 1
counts_atac1 = np.array(atac[:, regions].X.todense())
counts_atac1 = scmomat.preprocess(counts_atac1, modality = "ATAC")
counts_atacs = [counts_atac1, None]

# scRNA-seq of batch 1
counts_rna2 = np.array(rna[:, genes].X.todense())
counts_rna2 = scmomat.preprocess(counts_rna2, modality = "RNA", log = False)
counts_rnas = [None, counts_rna2,]

# CALCULATE THE PSEUDO-COUNT MATRIX
# Load the gene by region association matrix
# A = sp.load_npz(os.path.join(data_dir, 'GxR.npz'))
# A = np.array(A.todense())




# Calculate the pseudo-scRNA-seq matrix of batch 1
counts_rnas[0] = counts_atacs[0] @ A.T
# Binarize
counts_rnas[0] = (counts_rnas[0]!=0).astype(int)


# CREATE THE COUNTS OBJECT
counts = {"feats_name": feats_name, "nbatches": n_batches, "rna":counts_rnas, "atac": counts_atacs}


# ## Step 2: training scMoMaT
# The hyper-parameter includes `lamb` for weight lambda in the objective function, `K` for the number of latent dimensions. The default parameter (`lamb = 0.001`, `K = 30`) works for most of the real datasets.

# In[3]:


#------------------------------------------------------------------------------------------------------------------------------------
# NOTE: Number of latent dimensions, key hyper-parameter, 20~30 works for most of the cases.
K = 30
#------------------------------------------------------------------------------------------------------------------------------------
# NOTE: Here we list other parameters in the function for illustration purpose, most of these parameters are set as default value.
# weight on regularization term, default value
lamb = 0.001 
# number of total iterations, default value
T = 4000
# print the result after each ``interval'' iterations, default value
interval = 1000
# batch size for each iteraction, default value
batch_size = 0.1
# learning rate, default value
lr = 1e-2
# random seed, default value
seed = 0
# running device, can be CPU or GPU
device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")

#------------------------------------------------------------------------------------------------------------------------------------

start_time = time.time()
model = scmomat.scmomat_model(counts = counts, K = K, batch_size = batch_size, interval = interval, lr = lr, lamb = lamb, seed = seed, device = device)
losses = model.train_func(T = T)
end_time = time.time()
print("running time: " + str(end_time - start_time))

# Plot loss function
x = np.linspace(0, T, int(T/interval)+1)
plt.plot(x, losses)


# In[4]:


# # read in the ground truth labels
# labels = []
# labels_batches = []
# for batch in range(n_batches):
#     labels.append(pd.read_csv(os.path.join(data_dir, 'meta_c' + str(batch + 1) + '.csv'), index_col=0)["cluster (remapped)"].values.squeeze()[::subsample])        
#     labels_batches.append(np.array(["batch " + str(batch)] * len(labels[-1]), dtype = object))


# NOTE: Plot the result before post-processing
# umap_op = UMAP(n_components = 2, n_neighbors = 30, min_dist = 0.2, random_state = 0) 
zs = model.extract_cell_factors()
# x_umap = umap_op.fit_transform(np.concatenate(zs, axis = 0))

# scmomat.plot_latent(x_umap, annos = np.concatenate(labels), batches = np.concatenate(labels_batches), mode = "separate", save = None, figsize = (15,30), axis_label = "UMAP", markerscale = 6, s = 5, label_inplace = True)
# scmomat.plot_latent(x_umap, annos = np.concatenate(labels_batches), mode = "joint", save = None, figsize = (15,10), axis_label = "UMAP", markerscale = 6, s = 5, label_inplace = True)


inte1 = sc.AnnData(zs[1])
inte2 = sc.AnnData(zs[0])
inte1.obs_names = rna.obs_names
inte2.obs_names = atac.obs_names
inte1.obs['cell_type'] = rna.obs['cell_type']
inte2.obs['cell_type'] = atac.obs['cell_type']


path = '../results/'+ data_id
if not os.path.exists(path):
    os.makedirs(path)

inte1.write(os.path.join(path, "scmomat_rna.h5ad"), compression="gzip")
inte2.write(os.path.join(path, "scmomat_atac.h5ad"), compression="gzip")

