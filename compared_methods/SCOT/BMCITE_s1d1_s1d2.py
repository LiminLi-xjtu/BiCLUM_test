#!/usr/bin/env python
# coding: utf-8

import sys
import src.utils as ut
import src.evals as evals
import src.scot2 as sc
import numpy as np
import os
import anndata
import pandas as pd



# In[2]:

####################################################################
numk = 50
data_id = "BMCITE"
####################################################################

dataset_dir = '../../data/bmcite/'
dataset_name = 's1d1_s1d2'
dataset_type = 'combine' #'same_name'

RNA_data = anndata.read(os.path.join(dataset_dir, 'GASM/bm.rna.' + dataset_name +'.h5ad'))
ADT_data = anndata.read(os.path.join(dataset_dir, 'GASM/bm.adt.' + dataset_name +'.h5ad'))
# cite_graph = pd.read_csv(os.path.join(dataset_dir, 'GASM/' + 'graph_' + dataset_type + '.csv'))

data1 = RNA_data.obsm['X_pca'] # [:, cite_graph.values[:, 1]]
data2 = ADT_data.obsm['X_apca'] # [:, cite_graph.values[:, 2]]

X = data1
y = data2

# initialize SCOT object
scot=sc.SCOT(X, y)
# call the alignment with z-score normalization
X_aligned, y_aligned = scot.align( k=numk, e=5e-3,  normalize=True, norm="l2")

inte = []
inte.append(X_aligned)
inte.append(y_aligned[0])

scot_inte = dict({"inte": inte})


path = '../results/'+ data_id + '/' + dataset_name
if not os.path.exists(path):
    os.makedirs(path)

np.save(os.path.join(path, 'scot.npy'), scot_inte)
