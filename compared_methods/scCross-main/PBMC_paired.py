#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import anndata
import scanpy as sc
import sccross
import pandas as pd
from matplotlib import rcParams
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
import os

# # Read data

# In[ ]:

dataset_name = 'PBMC_paired'
rcParams["figure.figsize"] = (4, 4)
rna = anndata.read_h5ad("./data/" +dataset_name +"/rna_preprocessed.h5ad")
atac = anndata.read_h5ad("./data/" +dataset_name + "/atac_preprocessed.h5ad")

rna.obs['domain'] = 'scRNA-seq'
atac.obs['domain'] = 'scATAC-seq'
# # Configure data

# In[ ]:


sccross.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer = 'counts',
     use_rep="X_pca"
)

sccross.models.configure_dataset(
    atac, "NB", use_highly_variable=False,
    use_rep="X_lsi"
)


# # MNN prior

# In[ ]:


sccross.data.mnn_prior([rna,atac])


# # Training

# In[ ]:


cross = sccross.models.fit_SCCROSS(
    {"rna": rna, "atac": atac},
    fit_kws={"directory": "sccross"}
)





# # Integration benchmark

# In[ ]:


rna.obsm["X_cross"] = cross.encode_data("rna", rna)
atac.obsm["X_cross"] = cross.encode_data("atac", atac)


# In[ ]:

rna.obsm["X_cross"] = cross.encode_data("rna", rna)
atac.obsm["X_cross"] = cross.encode_data("atac", atac)

inte1 = sc.AnnData(rna.obsm['X_cross'])
inte2 = sc.AnnData(atac.obsm['X_cross'])
inte1.obs_names = rna.obs_names
inte2.obs_names = atac.obs_names
inte1.obs['cell_type'] = rna.obs['cell_type']
inte2.obs['cell_type'] = atac.obs['cell_type']


path = '../results/'+ dataset_name
if not os.path.exists(path):
    os.makedirs(path)

inte1.write(os.path.join(path, "sccross_rna.h5ad"), compression="gzip")
inte2.write(os.path.join(path, "sccross_atac.h5ad"), compression="gzip")


# # Save model

# In[ ]:


cross.save(os.path.join(path, "cross.dill"))