#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import random
random.seed(1)
import pandas as pd
from load_result import *
from evaluation import evaluate



dataset_name = "BMMC_s1d1"
dataset_type = 'RNA_ATAC'
GAM_name = 'Signac'
paired=True


path = "../"
dataset_dir = path + '../data/'
result_dir = path + "results/" + dataset_name
eva_dir = path + "eva/" + dataset_name
vis_dir = path + 'vis/' + dataset_name

if not os.path.exists(vis_dir):
    os.makedirs(vis_dir)
if not os.path.exists(eva_dir):
    os.makedirs(eva_dir)

# adata = anndata.read_h5ad(result_dir + 'adata.h5ad')
# adata_gasm = anndata.read_h5ad(result_dir + 'adata_gasm.h5ad')
# annotations = np.load(result_dir + "/annotations.npz")
# anno_rna = annotations['anno_rna']
# anno_other = annotations['anno_other']
# anno_rna_gasm = annotations['anno_rna_gasm']
# anno_other_gasm = annotations['anno_other_gasm']

methods_raw = ['GLUE', 'JointMDS', 'MMDMA', 'Pamona', 'scDART', 'SCOT', 'scTopoGAN', 'Seurat', 'UnionCom']
methods_gam = ['raw_data', 'bindSC', 'LIGER', 'MultiMAP', 'uniPort', 'BiCLUM']
adata, anno_rna, anno_other, clu_rna, clu_other = load_result_raw(dataset_name, dataset_dir, result_dir, methods_raw)
adata_gasm, anno_rna_gasm, anno_other_gasm, clu_rna_gasm, clu_other_gasm = load_result_gasm(dataset_name, dataset_dir, result_dir, GAM_name, methods_gam)

adata.write(result_dir+"/adata.h5ad", compression='gzip')
adata_gasm.write(result_dir+"/adata_gasm.h5ad", compression='gzip')
np.savez_compressed(result_dir+"/annotations.npz", anno_rna=anno_rna, anno_other=anno_other,
                    anno_rna_gasm=anno_rna_gasm, anno_other_gasm=anno_other_gasm)


eva_metrics_raw = evaluate(adata, anno_rna, anno_other, paired)
eva_metrics_gasm = evaluate(adata_gasm, anno_rna_gasm, anno_other_gasm, paired)
eva_metrics = pd.concat([eva_metrics_raw, eva_metrics_gasm], axis=1)
eva_metrics.to_csv(eva_dir + '/eva_metrics_' + GAM_name + '.csv', index=True)

