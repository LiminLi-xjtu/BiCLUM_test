#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import anndata
import MultiMAP
import numpy as np
import os
import pandas as pd

sc.settings.set_figure_params(dpi=80)

####################################################################
data_id = "BMCITE"
dataset_dir = '../../data/bmcite/'
dataset_name = 's1d1_s1d2'
dataset_type = 'combine' #'combine' #'same_name'
ADT_data = anndata.read(os.path.join(dataset_dir, 'GASM/bm.adt.' + dataset_name +'.h5ad'))
RNA_data = anndata.read(os.path.join(dataset_dir, 'GASM/bm.rna.' + dataset_name +'.h5ad'))
cite_graph = pd.read_csv(os.path.join(dataset_dir, 'GASM/' + 'graph_' + dataset_type + '.csv'))
####################################################################


ADT_data.obs['source'] = 'ADT'
RNA_data.obs['source'] = 'RNA'


# In[3]:


[RNA_data.shape, ADT_data.shape]


rna = RNA_data[:, cite_graph.values[:, 1]]
adt = ADT_data[:, cite_graph.values[:, 2]]
# In[5]:


adata = MultiMAP.Integration([rna, adt], ['X_pca', 'X_apca'])
inte_data = sc.AnnData(adata.obsm['X_multimap'])
inte_data.obs_names = adata.obs_names
inte_data.obs['batch'] = adata.obs['source']
inte_data.obs['cell_type'] = adata.obs['cell_type'].astype(str)


# inte = sc.AnnData(inte_data.obsm['X_multimap'])
# inte.obs_names = inte_data.obs_names
# inte.obs['batch'] = inte_data.obs['source']
# inte.obs['cell_type'] = inte_data.obs['cell_type'].astype(str)

path = '../results/'+ data_id + '/MultiMAP'
if not os.path.exists(path):
    os.makedirs(path)

inte_data.write_h5ad(os.path.join(path, dataset_name + '_' + dataset_type + ".h5ad"))

