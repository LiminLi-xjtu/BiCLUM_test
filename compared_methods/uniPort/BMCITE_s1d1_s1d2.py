#!/usr/bin/env python
# coding: utf-8


import uniport as up
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import os
import sys

def main():
    ####################################################################
    data_id = "bmcite"
    dataset_dir = '../../data/bmcite/'
    dataset_name = 's1d1_s1d2'

    adata_adt = anndata.read(os.path.join(dataset_dir, 'GASM/bm.adt.' + dataset_name + '.h5ad'))
    adata_rna = anndata.read(os.path.join(dataset_dir, 'GASM/bm.rna.' + dataset_name + '.h5ad'))
    cite_graph = pd.read_csv(os.path.join(dataset_dir, 'GASM/' + 'graph.csv'))
    ####################################################################




    adata_adt.obs['domain_id'] = 0
    adata_adt.obs['domain_id'] = adata_adt.obs['domain_id'].astype('category')
    adata_adt.obs['source'] = 'ADT'

    adata_rna.obs['domain_id'] = 1
    adata_rna.obs['domain_id'] = adata_rna.obs['domain_id'].astype('category')
    adata_rna.obs['source'] = 'RNA'

    # Concatenate scADT-seq and scRNA-seq with common genes using `AnnData.concatenate`.

    adata_cm = adata_adt.concatenate(adata_rna, join='inner', batch_key='domain_id')


    # Preprocess data with common genes. Select 2,000 highly variable common genes. \

    # sc.pp.highly_variable_genes(adata_cm, n_top_genes=2000, flavor="seurat_v3")
    # sc.pp.normalize_total(adata_cm)
    # sc.pp.log1p(adata_cm)
    # sc.pp.highly_variable_genes(adata_cm, n_top_genes=2000, inplace=False, subset=True)
    # up.batch_scale(adata_cm)
    # sc.pp.scale(adata_cm)
    print(adata_cm.obs)


    # Preprocess scRNA-seq data. Select 2,000 highly variable genes as RNA specific.

    # # sc.pp.highly_variable_genes(adata_rna, n_top_genes=2000, flavor="seurat_v3")
    # sc.pp.normalize_total(adata_rna)
    # sc.pp.log1p(adata_rna)
    # sc.pp.highly_variable_genes(adata_rna, n_top_genes=2000, inplace=False, subset=True)
    # up.batch_scale(adata_rna)
    # # sc.pp.scale(adata_rna)
    #
    #
    # # Preprocess scADT-seq data.
    # # Select 2,000 highly variable genes as ADT speicifc.
    #
    # # sc.pp.highly_variable_genes(adata_adt, n_top_genes=2000, flavor="seurat_v3")
    # sc.pp.normalize_total(adata_adt)
    # sc.pp.log1p(adata_adt)
    # sc.pp.highly_variable_genes(adata_adt, n_top_genes=2000, inplace=False, subset=True)
    # up.batch_scale(adata_adt)
    # sc.pp.scale(adata_adt)


    # Save the preprocessed data for integration.

    # ### Integration with specific genes and optimal transport loss
    # Integrate the scADT-seq and scRNA-seq data using both common and dataset-specific genes by `Run` function in uniport. The latent representations of data are stored in `adata.obs['latent']`.

    adata = up.Run(adatas=[adata_adt,adata_rna], adata_cm=adata_cm, lambda_s=1.0, iteration=10000, loss_type='MSE')
    adata.obs['cell_type'] = adata.obs['cell_type'].astype(str)

    path = '../results/'+ data_id + '/uniPort'
    if not os.path.exists(path):
        os.makedirs(path)

    adata.write_h5ad(os.path.join(path, dataset_name + '.h5ad'),compression='gzip')

main()
