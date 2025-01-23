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
data_id = "Kidney"
# GAM = input("please input GAM ID:")
GAM = 'ArchR'
dataset_dir = '../../data/kidney/'
atac_genes = anndata.read(os.path.join(dataset_dir, 'GASM/scGAM_' + GAM + '.h5ad'))
rna = anndata.read(os.path.join(dataset_dir, 'GASM/rna.h5ad'))
atac_peaks = anndata.read(os.path.join(dataset_dir, 'GASM/atac.h5ad'))
####################################################################
if GAM == 'ArchR':
    name_atac_genes_ = pd.Series(atac_genes.obs_names)
    name_atac_peaks_ = pd.Series(atac_peaks.obs_names)
    name_atac_peaks_dict = pd.Series(name_atac_peaks_.index, index=name_atac_peaks_.values).to_dict()
    match_positions = name_atac_genes_.map(name_atac_peaks_dict)
    atac_peaks = atac_peaks[match_positions, :]

atac_genes.obs['source'] = 'GAM'
rna.obs['source'] = 'RNA'
atac_peaks.obs['source'] = 'ATAC'

# As part of its operation, MultiMAP will compute PCA reductions of combinations of the input data. In contrast to standard operating procedures for RNA analysis, it is recommended to avoid highly variable gene filtering, only removing unexpressed/scarcely present genes. This grants the method as much overlapping information as possible between the datasets to try to find commonalities.
#
# There are two forms of ATAC data imported - a peak matrix, and a gene space. It is important to have a shared gene space between the processed datasets to enable shared PCA computation. In the case of ATAC data, converting peaks to genes can be performed with tools such as [SnapATAC](https://github.com/r3fang/SnapATAC) or [Signac](https://satijalab.org/signac/).

# In[3]:


[rna.shape, atac_peaks.shape, atac_genes.shape]


# The input to MultiMAP is an AnnData object per dataset to integrate, with a log-normalised gene space in `.X` and a primary dimensionality reduction somewhere in `.obsm`. Once again, using a large gene space is encouraged to capture as much information about the dataset as possible. The data in `.X` should not be scaled, as the output object will be created by concatenating the individual inputs. As such, this will make the result better for downstream analysis - scaled expression data is less informative.
#
# We're going to compute a standard PCA for the RNA data, and represent the ATAC peak space as LSI coordinates after a TF-IDF transformation. A helper function to perform this is included as `Multimap.TFIDF_LSI()`. Once we compute this dimensionality reduction, we'll copy it over to the gene space representation of the ATAC data. We'll also briefly use a helper object to compute a PCA for the RNA data, as PCA computation requires the data to be scaled, and copy the resulting coordinates to the object with un-scaled data.

# In[4]:


MultiMAP.TFIDF_LSI(atac_peaks)
atac_genes.obsm['X_lsi'] = atac_peaks.obsm['X_lsi'].copy()

rna_pca = rna.copy()
sc.pp.scale(rna_pca)
sc.pp.pca(rna_pca)
rna.obsm['X_pca'] = rna_pca.obsm['X_pca'].copy()


# Now that we have all our datasets present with a log-normalised gene space in `.X` and a primary dimensionality reduction in `.obsm`, we can run `MultiMAP.Integration()`. The function accepts a list of the AnnData objects on input, along with a second list with the names of the primary dimensionality reductions to use (as we're using `X_pca` for our RNA and `X_lsi` for our ATAC).

# In[5]:


adata = MultiMAP.Integration([rna, atac_genes], ['X_pca', 'X_lsi'])

inte = sc.AnnData(adata.obsm['X_multimap'])
inte.obs_names = adata.obs_names
inte.obs['batch'] = adata.obs['source']
inte.obs['cell_type'] = adata.obs['cell_type'].astype(str)

# The output is an AnnData object with the MultiMAP embedding stored in `.obsm['X_multimap']`, plus corresponding neighbourhood graphs in `.obsp`.



path = '../results/'+ data_id + '/MultiMAP'
if not os.path.exists(path):
    os.makedirs(path)

inte.write_h5ad(os.path.join(path, GAM + '.h5ad'), compression='gzip')


