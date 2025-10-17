import numpy as np
import os
from process_data import processing
from train_data import training
import scanpy as sc
import anndata

# conda install -c conda-forge -c bioconda scglue  # CPU only
# conda install -c conda-forge -c bioconda scglue pytorch-gpu  # With GPU support
# conda install -c pytorch faiss-cpu
# conda install -c "pytorch/label/nightly" faiss-gpu

species = 'Human' #'Human', "Mouse

data_id = "PBMC_paired"
dataset_dir ='../../data/PBMC_paired/'
rna = anndata.read(os.path.join(dataset_dir, 'rna.h5ad'))
atac = anndata.read(os.path.join(dataset_dir, 'atac.h5ad'))
rna.obs['batch'] = 'RNA_c'
atac.obs['batch'] = 'ATAC_c'




# stpe 1: processing data
processing(rna, atac, data_id, species)

# stpe 2: training data
rna_res, atac_res = training(data_id)

inte1 = sc.AnnData(rna_res.obsm['X_glue'])
inte2 = sc.AnnData(atac_res.obsm['X_glue'])
inte1.obs_names = rna_res.obs_names
inte2.obs_names = atac_res.obs_names
inte1.obs['cell_type'] = rna_res.obs['cell_type']
inte2.obs['cell_type'] = atac_res.obs['cell_type']
inte1.obs['batch'] = rna.obs['batch']
inte2.obs['batch'] = atac.obs['batch']


path = '../results/'+ data_id
if not os.path.exists(path):
    os.makedirs(path)

inte1.write(os.path.join(path, "glue_rna.h5ad"))
inte2.write(os.path.join(path, "glue_atac.h5ad"))



