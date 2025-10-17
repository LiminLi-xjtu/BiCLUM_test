#!/usr/bin/env python
# coding: utf-8


from pamona import Pamona
import numpy as np
import os
import anndata
import pandas as pd




data_id = "BMCITE_s1d1_s1d2"
dataset_dir = '../../data/BMCITE_s1d1_s1d2/'
dataset_name = 's1d1_s1d2'
dataset_type = 'combine' #'same_name'



RNA_data = anndata.read(os.path.join(dataset_dir, 'GASM/bm.rna.' + dataset_name +'.h5ad'))
ADT_data = anndata.read(os.path.join(dataset_dir, 'GASM/bm.adt.' + dataset_name +'.h5ad'))

data1 = RNA_data.obsm['X_pca'] # [:, cite_graph.values[:, 1]]
data2 = ADT_data.obsm['X_apca'] # [:, cite_graph.values[:, 2]]

X1 = data1
X2 = data2
data = [X1, X2]



Pa = Pamona.Pamona(n_neighbors=30, Lambda=1)

integrated_data, T = Pa.run_Pamona(data)


Pamona_inte = dict({"inte": integrated_data})


path = './results/'+ data_id + '/' + dataset_name
if not os.path.exists(path):
    os.makedirs(path)

np.save(os.path.join(path, 'Pamona.npy'), Pamona_inte)



