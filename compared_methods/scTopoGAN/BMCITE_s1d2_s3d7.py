import numpy as np
import pandas as pd
import anndata
import os
import pandas
import scanpy as sc

from scTopoGAN_Functions import get_TopoAE_Embeddings, run_scTopoGAN


data_id = "BMCITE"

dataset_dir = '../../data/bmcite/'
dataset_name = 's1d2_s3d7'
dataset_type = 'combine' #'same_name'

RNA_data = anndata.read(os.path.join(dataset_dir, 'GASM/bm.rna.' + dataset_name +'.h5ad'))
ADT_data = anndata.read(os.path.join(dataset_dir, 'GASM/bm.adt.' + dataset_name +'.h5ad'))
# cite_graph = pd.read_csv(os.path.join(dataset_dir, 'GASM/' + 'graph_' + dataset_type + '.csv'))

RNA = sc.AnnData(RNA_data.obsm['X_pca'].copy()) # [:, cite_graph.values[:, 1]]
ADT = sc.AnnData(ADT_data.obsm['X_apca'].copy())
RNA.obs_names = RNA_data.obs_names
ADT.obs_names = ADT_data.obs_names

# RNA = np.load('/Users/guoyin/Desktop/UZH/contrastive_learning/method/BC6/data/PBMC1/PBMC1_data2.npy', allow_pickle=True)
# ATAC = np.load('/Users/guoyin/Desktop/UZH/contrastive_learning/method/BC6/data/PBMC1/PBMC1_data3.npy', allow_pickle=True)


# Step 1: Get TopoAE embeddings
# set topology_regulariser_coefficient between 0.5 to 3.0

target_latent = get_TopoAE_Embeddings(Manifold_Data = RNA, batch_size=50, autoencoder_model="MLPAutoencoder_PBMC",
                                      AE_arch = [50, 32, 32, 8], topology_regulariser_coefficient=3, initial_LR=0.001)

source_latent = get_TopoAE_Embeddings(Manifold_Data = ADT, batch_size=50, autoencoder_model="MLPAutoencoder_PBMC",
                                      AE_arch = [50, 32, 32, 8], topology_regulariser_coefficient=0.5, initial_LR=0.001)

## Step 2: Manifold alignment using scTopoGAN
source_aligned = run_scTopoGAN(source_latent, target_latent, source_tech_name="ADT", target_tech_name="RNA",
                               batch_size=512, topology_batch_size=1000, total_epochs=1001, num_iterations=20,
                               checkpoint_epoch=100, g_learning_rate=1e-4, d_learning_rate=1e-4, path_prefix="Results")
# source_aligned = run_scTopoGAN(source_latent, target_latent, source_tech_name="ADT", target_tech_name="RNA",
#                                batch_size=512, topology_batch_size=1000, total_epochs=1, num_iterations=1,
#                                checkpoint_epoch=1, g_learning_rate=1e-4, d_learning_rate=1e-4, path_prefix="Results")

inte = []
inte.append(target_latent.to_numpy())
inte.append(source_aligned.to_numpy())

scTopoGAN_inte = dict({"inte": inte})

path = '../results/'+ data_id + '/' + dataset_name
if not os.path.exists(path):
    os.makedirs(path)

np.save(os.path.join(path, 'scTopoGAN.npy'), scTopoGAN_inte)

# source_aligned.to_csv('ATAC_aligned.csv')