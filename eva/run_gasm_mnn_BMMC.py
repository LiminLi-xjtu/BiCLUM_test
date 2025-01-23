import scanpy as sc
import numpy as np
import os
import matplotlib.pyplot as plt
from datasets import load_data
from sklearn.decomposition import PCA
import pandas as pd
import os
import plotly.graph_objects as go
from sklearn.cluster import KMeans
from scipy.stats import pearsonr

from config import load_config
import argparse
import torch
import anndata as ad
from calculate_mnn import *

config_args = load_config('./config/BMMC_s1d1')
args = argparse.Namespace(**config_args)

dataset_name = "BMMC_s1d1"
dataset_type = 'RNA_ATAC'
args.matched = True

dataset_name = "BMCITE_s1d1_s1d2"
dataset_batch = 's1d1_s1d2'
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

k_mnn = [10, 30, 50, 80, 100, 150, 200, 300, 400, 500]
k_mnn_str = [str(k) for k in k_mnn]
GAM_name = ['ArchR', 'Cicero_unnorm', 'cisTopic', 'dsc', 'MAESTRO','Signac', 'snapATAC2'] #'ArchR', 'Cicero_unnorm', 'cisTopic', 'dsc', 'MAESTRO',


col_names = ['mnn_pairs', 'correct_ratio', 'mnn_corr']
df = pd.DataFrame(index=k_mnn_str, columns=col_names)


path = "../"
dataset_dir = path + '../data/processed_data/'
result_dir = path + "results/" + dataset_name
eva_dir = path + "eva_mix/" + dataset_name
vis_dir = path + 'vis/' + dataset_name
eva_pearson_cell_dir = path + "eva_mix/" + dataset_name + '/eva_pearson_cell/'
eva_pearson_gene_dir = path + "eva_mix/" + dataset_name + '/eva_pearson_gene/'

if not os.path.exists(vis_dir):
    os.makedirs(vis_dir)
if not os.path.exists(eva_dir):
    os.makedirs(eva_dir)
if not os.path.exists(eva_pearson_cell_dir):
    os.makedirs(eva_pearson_cell_dir)
if not os.path.exists(eva_pearson_gene_dir):
    os.makedirs(eva_pearson_gene_dir)


use_cuda = torch.cuda.is_available()
device = torch.device('cuda:0' if use_cuda else 'cpu')
args.device = device




for i in range(len(GAM_name)):

    args.GAM_name = GAM_name[i]
    print(f"load data: {args.GAM_name}")
    gam, rna, atac, label, meta_batch, meta_omic = load_data(args.dataset_name, args.dataset_dir, args.GAM_name)


    GAM_data, RNA_data, args = data_input(args, meta_batch)

    gam_names = GAM_data.obs_names
    rna_names = RNA_data.obs_names

    gam_hvg = GAM_data.var.query("highly_variable").index
    rna_hvg = RNA_data.var.query("highly_variable").index
    data_mnn = [GAM_data.X, RNA_data.X]



    print('correlation between cells')
    pearson_cell = np.zeros((np.shape(GAM_data)[0], 1))
    for p in range(np.shape(GAM_data)[0]):
        # pearson_cell[p, 0] = np.corrcoef(GAM_data.X[p, :], RNA_data_.X[p, :])[0, 1]
        # pearson_cell[p, 0], _ = pearsonr(np.array(GAM_data.X[p, :]), np.array(RNA_data_.X[p, :]))
        pearson_cell[p, 0], _ = pearsonr(gam_x[p, :], rna_x[p, :])
    pearson_cell_avg = np.mean(pearson_cell)

    print('correlation between gene')
    pearson_gene = np.zeros((np.shape(GAM_data)[1], 1))
    for q in range(np.shape(GAM_data)[1]):
        # pearson_gene[q, 0] = np.corrcoef(GAM_data.X[:, q], RNA_data_.X[:, q])[0, 1]
        pearson_gene[q, 0], _ = pearsonr(gam_x[:, q], rna_x[:, q])
    pearson_gene_avg = np.mean(pearson_gene)

    df_pearson_cell = pd.DataFrame(pearson_cell, columns=['Corr'])
    df_pearson_gene = pd.DataFrame(pearson_gene, columns=['Corr'])

    df_pearson_cell.to_csv(eva_pearson_cell_dir + '_3.csv', index=True)
    df_pearson_gene.to_csv(eva_pearson_gene_dir + '_3.csv',
                            index=True)

    for j in range(len(k_mnn)):

        print(f"construct positive pairs:", {k_mnn[j]})

        args.neighbors_mnn = k_mnn[j]
        if args.dataset_type == 'RNA_ATAC':
            anchor_mnn_raw = construct_pairs_mnn_cell([GAM_data.X, RNA_data.X], args, meta_batch)
        else:
            anchor_mnn_raw = construct_pairs_mnn_cell([GAM_data[:,gam_names], RNA_data[:,rna_names]], args, meta_batch)

        num_anchor, mnn_num_type, mnn_ratio_type = cal_mnn(anchor_mnn_raw, label, args)

        mnn_pearson_vec = np.zeros((1, num_anchor))

        for p in range(num_anchor):
            mnn_pearson_vec[0, p] = np.corrcoef(data_mnn[0][anchor_mnn_raw[p, 0]], data_mnn[1][anchor_mnn_raw[p, 1]])[
                0, 1]
        mnn_pearson = mnn_pearson_vec.mean()



        df.iloc[j, 0] = mnn_num_type
        df.iloc[j, 1] = mnn_ratio_type
        df.iloc[j, 2] = mnn_pearson



    del GAM_data, RNA_data, data_mnn, anchor_mnn_raw


    df.to_csv(eva_dir + '/eva_mnn_' + args.dataset_name + '_' +  args.GAM_name + '.csv', index=True)
    # mnn_pearson_vec.to_csv(eva_dir + '/eva_mnn_corr' + args.dataset_name + '_' +  args.GAM_name + '.csv', index=True)
        

