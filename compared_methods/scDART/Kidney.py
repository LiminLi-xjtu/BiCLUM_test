#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys

sys.path.append('../')

import numpy as np
import pandas as pd

import torch
from sklearn.decomposition import PCA

import scDART.utils as utils
import scDART.TI as ti
import scDART
import anndata
import os

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

####################################################################

data_id = "Kidney"
path = "../"
dataset_dir = path + '../data/'
adata_rna = anndata.read(os.path.join(dataset_dir, data_id+'/GASM/rna.h5ad'))
adata_atac = anndata.read(os.path.join(dataset_dir, data_id+'/GASM/atac.h5ad'))
reg = anndata.read('../datasets/scDART/data/' + data_id + '.h5ad')
####################################################################

# all in one
seeds = [0]
latent_dim = 32
learning_rate = 3e-4
n_epochs = 500
use_anchor = False
reg_d = 1
reg_g = 1
reg_mmd = 1
ts = [30, 50, 70]
use_potential = True



gene_names = reg.obs_names
peak_names = reg.var_names
coarse_reg = reg.T.X#.todense()

counts_atac = adata_atac[:, peak_names].X#.todense()
counts_rna = adata_rna[:, gene_names].X#.todense()


scDART_op = scDART.scDART(n_epochs=n_epochs, latent_dim=latent_dim, \
                          ts=ts, use_anchor=use_anchor, use_potential=use_potential, k=10, \
                          reg_d=1, reg_g=1, reg_mmd=1, l_dist_type='kl', seed=seeds[0], \
                          device=torch.device('cuda' if torch.cuda.is_available() else 'cpu'))

scDART_op = scDART_op.fit(rna_count=counts_rna, atac_count=counts_atac, reg=coarse_reg, rna_anchor=None,
                          atac_anchor=None)
z_rna, z_atac = scDART_op.transform(rna_count=counts_rna, atac_count=counts_atac, rna_anchor=None,
                                    atac_anchor=None)

inte = []
inte.append(z_rna)
inte.append(z_atac)

scDART_inte = dict({"inte": inte})


path = '../results/'+ data_id
if not os.path.exists(path):
    os.makedirs(path)

np.save(os.path.join(path, 'scDART.npy'), scDART_inte)

