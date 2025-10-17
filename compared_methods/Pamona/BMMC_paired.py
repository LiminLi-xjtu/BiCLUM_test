#!/usr/bin/env python
# coding: utf-8


from pamona import Pamona
import numpy as np
import os
import anndata


data_id = "BMMC_paired"

path = '../datasets/' + data_id
data1 = anndata.read(os.path.join(path, "raw_data_rna.h5ad"))
data2 = anndata.read(os.path.join(path, "raw_data_atac.h5ad"))
X1 = data1.X
X2 = data2.X
data = [X1, X2]



Pa = Pamona.Pamona(n_neighbors=30, Lambda=1)

integrated_data, T = Pa.run_Pamona(data)


Pamona_inte = dict({"inte": integrated_data})


path = './results/'+ data_id
if not os.path.exists(path):
    os.makedirs(path)

np.save(os.path.join(path, 'Pamona.npy'), Pamona_inte)



