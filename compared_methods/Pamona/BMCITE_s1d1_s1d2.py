#!/usr/bin/env python
# coding: utf-8


from pamona import Pamona
import numpy as np
import os
import anndata
import pandas as pd




data_id = "BMCITE"
dataset_dir = '../../data/bmcite/'
dataset_name = 's1d1_s1d2'
dataset_type = 'combine'

# data1 = anndata.read(os.path.join(path, "raw_data_rna.h5ad"))
# data2 = anndata.read(os.path.join(path, "raw_data_atac.h5ad"))

RNA_data = anndata.read(os.path.join(dataset_dir, 'GASM/bm.rna.' + dataset_name +'.h5ad'))
ADT_data = anndata.read(os.path.join(dataset_dir, 'GASM/bm.adt.' + dataset_name +'.h5ad'))
# cite_graph = pd.read_csv(os.path.join(dataset_dir, 'GASM/' + 'graph_' + dataset_type + '.csv'))

data1 = RNA_data.obsm['X_pca'] # [:, cite_graph.values[:, 1]]
data2 = ADT_data.obsm['X_apca'] # [:, cite_graph.values[:, 2]]

X1 = data1
X2 = data2
data = [X1, X2]

# path = '../datasets/' + data_id
# data = np.load(os.path.join(path, 'rawdata.npy'), allow_pickle=True).item()
# data1 = data['exp'][0]
# data2 = data['exp'][1]
# data = [data1, data2]


Pa = Pamona.Pamona(n_neighbors=50, Lambda=10)

integrated_data, T = Pa.run_Pamona(data)


Pamona_inte = dict({"inte": integrated_data})


path = '../results/'+ data_id + '/' + dataset_name
if not os.path.exists(path):
    os.makedirs(path)

np.save(os.path.join(path, 'Pamona.npy'), Pamona_inte)



