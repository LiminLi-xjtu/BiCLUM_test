#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import os
import random
random.seed(1)
import pandas as pd
from load_result import *
from evaluation import evaluate

dataset_name = "BMCITE_s1d2_s3d7"
dataset_batch = 's1d2_s3d7'
dataset_type = 'RNA_ADT'
paired = True

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
# annotations = np.load(result_dir + "/annotations.npz")
# anno_rna = annotations['anno_rna']
# anno_other = annotations['anno_other']

methods = ['raw_data', 'Pamona', 'SCOT', 'scTopoGAN', 'UnionCom'] + ['bindSC', 'LIGER', 'MultiMAP', 'Seurat', 'uniPort', 'BiCLUM']
adata, anno_rna, anno_other, clu_rna, clu_other = load_cite(dataset_name, dataset_dir, result_dir, methods, dataset_batch)
adata.write(result_dir+"/adata.h5ad", compression='gzip')
np.savez_compressed(result_dir+"/annotations.npz", anno_rna=anno_rna, anno_other=anno_other, clu_rna=clu_rna, clu_other=clu_other)

eva_metrics = evaluate(adata, anno_rna, anno_other, paired)
eva_metrics.to_csv(eva_dir + '/eva_metrics.csv', index=True)

