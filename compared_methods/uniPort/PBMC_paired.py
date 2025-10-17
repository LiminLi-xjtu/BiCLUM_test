#!/usr/bin/env python
# coding: utf-8


import uniport as up
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import os
import sys

####################################################################
data_id = "PBMC_paired"
# GAM = 'ArchR' #input("please input GAM ID:")
GAM = sys.argv[1]
dataset_dir = '../../data/PBMC_paired/' # '../../../../data/guoyin/1-comp_methods/datasets/8-PBMC1/'
adata_atac = anndata.read(os.path.join(dataset_dir, 'GASM/scGAM_' + GAM + '.h5ad'))
adata_rna = anndata.read(os.path.join(dataset_dir, 'GASM/rna.h5ad'))
####################################################################

adata_atac.obs['domain_id'] = 0
adata_atac.obs['domain_id'] = adata_atac.obs['domain_id'].astype('category')
adata_atac.obs['source'] = 'ATAC'

adata_rna.obs['domain_id'] = 1
adata_rna.obs['domain_id'] = adata_rna.obs['domain_id'].astype('category')
adata_rna.obs['source'] = 'RNA'

# Concatenate scATAC-seq and scRNA-seq with common genes using `AnnData.concatenate`.

adata_cm = adata_atac.concatenate(adata_rna, join='inner', batch_key='domain_id')


# Preprocess data with common genes. Select 2,000 highly variable common genes. \

# sc.pp.highly_variable_genes(adata_cm, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(adata_cm)
sc.pp.log1p(adata_cm)
sc.pp.highly_variable_genes(adata_cm, n_top_genes=2000, inplace=False, subset=True)
up.batch_scale(adata_cm)
# sc.pp.scale(adata_cm)
print(adata_cm.obs)


# Preprocess scRNA-seq data. Select 2,000 highly variable genes as RNA specific.

# sc.pp.highly_variable_genes(adata_rna, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(adata_rna)
sc.pp.log1p(adata_rna)
sc.pp.highly_variable_genes(adata_rna, n_top_genes=2000, inplace=False, subset=True)
up.batch_scale(adata_rna)
# sc.pp.scale(adata_rna)


# Preprocess scATAC-seq data.
# Select 2,000 highly variable genes as ATAC speicifc.

# sc.pp.highly_variable_genes(adata_atac, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(adata_atac)
sc.pp.log1p(adata_atac)
sc.pp.highly_variable_genes(adata_atac, n_top_genes=2000, inplace=False, subset=True)
up.batch_scale(adata_atac)
# sc.pp.scale(adata_atac)


# Save the preprocessed data for integration.

# ### Integration with specific genes and optimal transport loss
# Integrate the scATAC-seq and scRNA-seq data using both common and dataset-specific genes by `Run` function in uniport. The latent representations of data are stored in `adata.obs['latent']`.

adata = up.Run(adatas=[adata_atac,adata_rna], adata_cm=adata_cm, lambda_s=1.0, iteration=10000)
# adata.obs['cell_type'] = adata.obs['cell_type'].astype(str)
# adata.obs['orig.ident'] = adata.obs['orig.ident'].astype(str)


inte = sc.AnnData(adata.obsm['latent'])
inte.obs_names = adata.obs_names
inte.obs['batch'] = adata.obs['source'].astype(str)
inte.obs['cell_type'] = adata.obs['cell_type'].astype(str)
inte.obs['orig.ident'] = adata.obs['orig.ident'].astype(str)

path = '../results/'+ data_id + '/uniPort'
if not os.path.exists(path):
    os.makedirs(path)

inte.write_h5ad(os.path.join(path, GAM + '.h5ad'),compression='gzip')


