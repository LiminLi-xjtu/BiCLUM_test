# In[]
import sys, os

import numpy as np
from umap import UMAP
import time
import torch
import matplotlib.pyplot as plt
import pandas as pd  
import scipy.sparse as sp
from sklearn.decomposition import PCA
import scanpy as sc
import numpy as np

import scmomat.model as model
import scmomat.utils as utils
import scmomat.bmk as bmk
import scmomat.umap_batch as umap_batch

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
def lsi(counts, n_components = 30):
    from sklearn.feature_extraction.text import TfidfTransformer
    from sklearn.decomposition import TruncatedSVD

    tfidf = TfidfTransformer(norm='l2', sublinear_tf=True)
    normed_count = tfidf.fit_transform(counts)

    # perform SVD on the sparse matrix
    lsi = TruncatedSVD(n_components=n_components + 1, random_state=42)
    lsi_r = lsi.fit_transform(normed_count)

    lsi.explained_variance_ratio_

    X_lsi = lsi_r[:, 1:]
    return X_lsi

# In[]
# ------------------------------------------------------------------------------------------------------------------------------------------------------
#
#   NOTE: 1. Load dataset and running scmomat (without retraining, retraining see the third section)
#
# ------------------------------------------------------------------------------------------------------------------------------------------------------
# NOTE: read in dataset

n_batches = 2
counts_rnas = []
counts_proteins = []

################TODO:rna,adt,network###############

data_id = 'BMCITE_s1d2_s3d7'
species = 'Human' 
adata_rna = sc.read_h5ad('../../data/BMCITE_s1d2_s3d7/bm.rna.s1d2_s3d7_new.h5ad')
adata_adt = sc.read_h5ad('../../data/BMCITE_s1d2_s3d7/bm.adt.s1d2_s3d7_new.h5ad')
fea_pairs = pd.read_csv('../../data/BMCITE_s1d2_s3d7/graph.csv')

counts_rna = np.array(adata_rna.X.todense())
counts_rna = utils.preprocess(counts_rna, modality = "RNA", log = False)
counts_protein = np.array(adata_adt.X.todense())
counts_protein = utils.preprocess(counts_protein, modality = "RNA", log = True)

A = np.zeros((adata_rna.shape[1], adata_adt.shape[1]))
proteins = fea_pairs.values[:, 2]
genes = fea_pairs.values[:, 1]

# 建立 b 中元素到索引的映射
pos_map = {val: idx for idx, val in enumerate(np.array(adata_rna.var_names))}
# 找到 a 每个元素在 b 中的位置
position_genes = [pos_map[x] for x in genes]

# 建立 b 中元素到索引的映射
pos_map2 = {val: idx for idx, val in enumerate(np.array(adata_adt.var_names))}
# 找到 a 每个元素在 b 中的位置
position_proteins = [pos_map2[x] for x in proteins]

A[position_genes, position_proteins] = 1

counts_rnas = [counts_rna, None]
counts_proteins = [None, counts_protein]

# for batch in range(n_batches):

#     try:
#         # counts_adt = sp.load_npz(os.path.join(dir, 'RxC' + str(batch + 1) + ".npz"))
#         # mmwrite(os.path.join(dir, 'RxC' + str(batch + 1) + ".mtx"), counts_adt)
#         counts_adt = sp.load_npz(os.path.join(dir, 'RxC' + str(batch + 1) + ".npz")).toarray().T
#         counts_adt = utils.preprocess(counts_adt, modality = "adt")
#     except:
#         counts_adt = None
        
#     try:
#         # counts_rna = sp.load_npz(os.path.join(dir, 'GxC' + str(batch + 1) + ".npz"))
#         # mmwrite(os.path.join(dir, 'GxC' + str(batch + 1) + ".mtx"), counts_rna)
#         counts_rna = sp.load_npz(os.path.join(dir, 'GxC' + str(batch + 1) + ".npz")).toarray().T
#         counts_rna = utils.preprocess(counts_rna, modality = "RNA", log = False)
#     except:
#         counts_rna = None
    
#     try:
#         # counts_protein = sp.load_npz(os.path.join(dir, 'PxC' + str(batch + 1) + ".npz"))
#         # mmwrite(os.path.join(dir, 'PxC' + str(batch + 1) + ".mtx"), counts_protein)
#         # the log transform produce better results for the protein
#         counts_protein = sp.load_npz(os.path.join(dir, 'PxC' + str(batch + 1) + ".npz")).toarray().T
#         counts_protein = utils.preprocess(counts_protein, modality = "RNA", log = True)
#     except:
#         counts_protein = None
    
#     # preprocess the count matrix
#     counts_rnas.append(counts_rna)
#     counts_adts.append(counts_adt)
#     counts_proteins.append(counts_protein)

counts = {"rna":counts_rnas, "protein": counts_proteins}

# A1 = sp.load_npz(os.path.join(dir, 'GxP.npz')).toarray()
# A2 = sp.load_npz(os.path.join(dir, 'GxR.npz')).toarray()

# obtain the feature name
# genes = pd.read_csv(dir + "genes.txt", header = None).values.squeeze()
# regions = pd.read_csv(dir + "regions.txt", header = None).values.squeeze()
# proteins = pd.read_csv(dir + "proteins.txt", header = None).values.squeeze()

feats_name = {"rna": np.array(adata_rna.var_names), "protein": np.array(adata_adt.var_names)}
counts["feats_name"] = feats_name

counts["nbatches"] = n_batches
# In[]
# NOTE: Running scmomat
# weight on regularization term
lamb = 0.001
batchsize = 0.1
# running seed
seed = 0
# number of latent dimensions
K = 30
interval = 1000
T = 4000
lr = 1e-2

start_time = time.time()
model1 = model.scmomat_model(counts = counts, K = K, batch_size = batchsize, interval = interval, lr = lr, lamb = lamb, seed = seed, device = device)
losses1 = model1.train_func(T = T)
end_time = time.time()
print("running time: " + str(end_time - start_time))

# x = np.linspace(0, T, int(T/interval)+1)
# plt.plot(x, losses1)
# plt.yscale("log")

# torch.save(model1, result_dir + f'CFRM_{K}_{T}.pt')
# model1 = torch.load(result_dir + f'CFRM_{K}_{T}.pt')

# In[] Sanity check, the scales should be positive, A_assos should also be positive
for mod in model1.A_assos.keys():
    if mod != "shared":
        print(torch.min(model1.A_assos["shared"] + model1.A_assos[mod]).item())

for mod in model1.A_assos.keys():
    if mod != "shared":
        print(torch.mean(model1.A_assos["shared"] + model1.A_assos[mod]).item())

for mod in model1.A_assos.keys():
    if mod != "shared":
        print(torch.max(model1.A_assos["shared"] + model1.A_assos[mod]).item())

print(model1.scales)

# In[]
# NOTE: Plot the result before post-processing (no post-processing for pbmc)
umap_op = UMAP(n_components = 2, n_neighbors = 30, min_dist = 0.2, random_state = 0) 
zs = []
for batch in range(n_batches):
    z = model1.softmax(model1.C_cells[str(batch)].cpu().detach()).numpy()
    zs.append(z)



inte1 = sc.AnnData(zs[0])
inte2 = sc.AnnData(zs[1])
inte1.obs_names = adata_rna.obs_names
inte2.obs_names = adata_adt.obs_names
inte1.obs['cell_type'] = adata_rna.obs['cell_type']
inte2.obs['cell_type'] = adata_adt.obs['cell_type']


path = '../results/'+ data_id
if not os.path.exists(path):
    os.makedirs(path)

inte1.write(os.path.join(path, "scmomat_rna.h5ad"), compression="gzip")
inte2.write(os.path.join(path, "scmomat_adt.h5ad"), compression="gzip")