#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import anndata
import MultiMAP
import numpy as np
import os

sc.settings.set_figure_params(dpi=80)

####################################################################
data_id = "BMMC_paired"
GAM = 'ArchR'

path = "../"
dataset_dir = path + '../data/'
result_dir = path + "results/" + data_id


GAM_data = anndata.read(os.path.join(dataset_dir, data_id+'/GASM/scGAM_' + GAM + '.h5ad'))
rna = anndata.read(os.path.join(dataset_dir, data_id+'/GASM/rna.h5ad'))
ATAC_data = anndata.read(os.path.join(dataset_dir, data_id+'/GASM/atac.h5ad'))
####################################################################


GAM_data.obs['source'] = 'GAM'
rna.obs['source'] = 'RNA'
ATAC_data.obs['source'] = 'ATAC'

# As part of its operation, MultiMAP will compute PCA reductions of combinations of the input data. In contrast to standard operating procedures for RNA analysis, it is recommended to avoid highly variable gene filtering, only removing unexpressed/scarcely present genes. This grants the method as much overlapping information as possible between the datasets to try to find commonalities.
#
# There are two forms of ATAC data imported - a peak matrix, and a gene space. It is important to have a shared gene space between the processed datasets to enable shared PCA computation. In the case of ATAC data, converting peaks to genes can be performed with tools such as [SnapATAC](https://github.com/r3fang/SnapATAC) or [Signac](https://satijalab.org/signac/).

# In[3]:


[rna.shape, ATAC_data.shape, GAM_data.shape]

var_names_1 = set(ATAC_data.obs_names)
var_names_2 = set(GAM_data.obs_names)
common_var_names = var_names_1.intersection(var_names_2)
atac_peaks = ATAC_data[list(common_var_names),:]
atac_genes = GAM_data[list(common_var_names),:]


# adata_fm = GAM_data.concatenate(ATAC_data, join='outer', batch_key='source')
# atac_peaks = adata_fm[adata_fm.obs["source"] == 'GAM']
# atac_genes = adata_fm[adata_fm.obs["source"] == 'ATAC']


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

# GAM_name = 'ArchR'
# GAM_name = 'cisTopic'
# GAM_name = 'Signac'
# GAM_name = 'snapATAC2'
#
# adata = anndata.read_h5ad(os.path.join(GAM_name + ".h5ad"))
# inte = sc.AnnData(adata.obsm['X_multimap'])
# inte.obs_names = adata.obs_names
# inte.obs['batch'] = adata.obs['source']
# inte.obs['cell_type'] = adata.obs['cell_type'].astype(str)
# inte.write_h5ad(os.path.join('inte/' + GAM_name + ".h5ad"))

# The output is an AnnData object with the MultiMAP embedding stored in `.obsm['X_multimap']`, plus corresponding neighbourhood graphs in `.obsp`.

# In[6]:

# inte_all = adata.obsm['X_multimap']
# inte = []
# batch = np.unique(adata.obs['source'])
# for i in batch:
#     inte.append(inte_all[adata.obs['source']==i,:])

path = '../results/'+ data_id + '/MultiMAP'
if not os.path.exists(path):
    os.makedirs(path)


inte.write_h5ad(os.path.join(path, GAM + '_new.h5ad'))

# sc.pl.embedding(adata, 'X_multimap', color=['source','cell_type'])

