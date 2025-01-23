
import numpy as np
import anndata
from sklearn.preprocessing import LabelEncoder
import os

def load_clusters(anno_other, anno_rna):

    all_cell_types = list(set(anno_other.tolist() + anno_rna.tolist()))
    label_encoder = LabelEncoder()
    label_encoder.fit_transform(all_cell_types)
    clu_other = label_encoder.transform(anno_other)
    clu_rna = label_encoder.transform(anno_rna)

    return clu_rna, clu_other



def load_result_raw(dataset_name, dataset_dir, result_dir, methods_raw):

    RNA_data = anndata.read(os.path.join(dataset_dir, dataset_name + '/GASM/rna.h5ad'))
    Other_data = anndata.read(os.path.join(dataset_dir, dataset_name + '/GASM/atac.h5ad'))
    Other_data.obs['omic_id'] = 'Other'
    RNA_data.obs['omic_id'] = 'RNA'
    anno_rna, anno_other = RNA_data.obs['cell_type'], Other_data.obs['cell_type']
    clu_rna, clu_other = load_clusters(anno_other, anno_rna)

    adata = anndata.read_h5ad(os.path.join(result_dir, "Seurat.h5ad"))
    adata.obsm = {}
    adata.obs['cell_type'] = np.concatenate((anno_rna, anno_other), axis = 0)
    adata.obs['cluster'] = np.concatenate((clu_rna, clu_other), axis=0)
    adata.obs['omic_id'] = np.concatenate((RNA_data.obs['omic_id'], Other_data.obs['omic_id']), axis=0)

    for i in range(len(methods_raw)):

        method = methods_raw[i]
        print(method)

        if method == 'GLUE':
            res_Other = anndata.read_h5ad(os.path.join(result_dir, 'glue_atac.h5ad'))
            res_RNA = anndata.read_h5ad(os.path.join(result_dir, 'glue_rna.h5ad'))
            adata.obsm[method] = np.vstack((res_RNA.X, res_Other.X))

        elif method == 'Seurat':
            result = anndata.read_h5ad(os.path.join(result_dir, "Seurat.h5ad"))
            if 'X_pca' in result.obsm:
                emb = np.array(result.obsm['X_pca'])
            else:
                try:
                    emb = np.array(result.X)
                except:
                    emb = np.array(result.X.todense())
            adata.obsm[method] = emb
        else:
            result = np.load(os.path.join(result_dir, method + ".npy"), allow_pickle=True).item()
            adata.obsm[method] = np.vstack((result['inte'][0], result['inte'][1]))  # rna, other

    return adata, anno_rna, anno_other, clu_rna, clu_other



def load_result_gasm(dataset_name, dataset_dir, result_dir, GAM_name, methods_gam):

    RNA_data = anndata.read(os.path.join(dataset_dir, 'processed_data/' + dataset_name + '_' + GAM_name + "-rna-pp.h5ad"))
    Other_data = anndata.read(os.path.join(dataset_dir, 'processed_data/' + dataset_name + '_' + GAM_name + "-gam-pp.h5ad"))

    anno_other = Other_data.obs['cell_type']
    anno_rna = RNA_data.obs['cell_type']
    clu_rna, clu_other = load_clusters(anno_other, anno_rna)
    num_other, num_rna = len(clu_other), len(clu_rna)
    Other_data.obs['omic_id'] = 'Other'
    RNA_data.obs['omic_id'] = 'RNA'

    adata = anndata.concat([RNA_data, Other_data], join='outer')
    adata.obsm = {}
    adata.obs['cell_type'] = np.concatenate((anno_rna, anno_other), axis = 0)
    adata.obs['cluster'] = np.concatenate((clu_rna, clu_other), axis=0)
    adata.obs['omic_id'] = np.concatenate((RNA_data.obs['omic_id'], Other_data.obs['omic_id']), axis=0)


    for i in range(len(methods_gam)):

        method = methods_gam[i]

        if method == 'raw_data':
            adata.obsm[method] = np.vstack((RNA_data.obsm['X_pca'], Other_data.obsm['X_pca']))

        elif method == 'bindSC':
            result = anndata.read_h5ad(os.path.join(result_dir, "bindSC/" + GAM_name + ".h5ad"))
            emb = np.array(result.X.todense())
            adata.obsm[method] = emb  # [emb[0:num_rna, :], emb[num_rna:, :]]

        elif method == 'LIGER':
            result = anndata.read_h5ad(os.path.join(result_dir, "LIGER/" + GAM_name + ".h5ad"))
            emb = np.array(result.X.todense())
            adata.obsm[method] = np.vstack((emb[num_other:, :], emb[0:num_other, :]))

        elif method == 'MultiMAP':
            result = anndata.read_h5ad(os.path.join(result_dir, "MultiMAP/" + GAM_name + ".h5ad"))
            emb = result.X
            adata.obsm[method] = emb

        elif method == 'uniPort':
            result = anndata.read_h5ad(os.path.join(result_dir, "uniPort/" + GAM_name + ".h5ad"))
            emb = np.array(result.X) # .todense()
            adata.obsm[method] = np.vstack((emb[num_other:, :], emb[0:num_other, :]))

        elif method == 'BiCLUM':
            result = np.load(os.path.join(result_dir, "BiCLUM/" + GAM_name + '.npy'), allow_pickle=True).item()
            adata.obsm[method] = np.vstack((result['inte'][1], result['inte'][0]))

    return adata, anno_rna, anno_other, clu_rna, clu_other



def load_cite(dataset_name, dataset_dir, result_dir, methods, batch):

    RNA_data = anndata.read(os.path.join(dataset_dir, dataset_name + '/bm.rna.' + batch + '.h5ad'))
    Other_data = anndata.read(os.path.join(dataset_dir, dataset_name + '/bm.adt.' + batch + '.h5ad'))
    Other_data.obs['omic_id'] = 'Other'
    RNA_data.obs['omic_id'] = 'RNA'
    anno_rna, anno_other = RNA_data.obs['cell_type2'], Other_data.obs['cell_type2']
    clu_rna, clu_other = load_clusters(anno_other, anno_rna)
    num_other, num_rna = len(clu_other), len(clu_rna)

    adata = anndata.concat([RNA_data, Other_data], join='outer')
    adata.obs_names_make_unique()
    adata.obsm = {}
    adata.obs['cell_type'] = np.concatenate((anno_rna, anno_other), axis = 0)
    adata.obs['cluster'] = np.concatenate((clu_rna, clu_other), axis=0)
    adata.obs['omic_id'] = np.concatenate((RNA_data.obs['omic_id'], Other_data.obs['omic_id']), axis=0)

    for i in range(len(methods)):
        method = methods[i]
        print(method)

        if method == 'raw_data':
                adata.obsm[method] = np.vstack((RNA_data.obsm['X_pca'], Other_data.obsm['X_apca']))

        elif method == 'bindSC':
            result = anndata.read_h5ad(os.path.join(result_dir, "bindSC.h5ad"))
            adata.obsm[method] = np.array(result.X.todense())

        elif method == 'LIGER':
            result = anndata.read_h5ad(os.path.join(result_dir, "LIGER.h5ad"))
            emb = np.array(result.X.todense())
            adata.obsm[method] = np.vstack((emb[num_other:, :], emb[0:num_other, :]))

        elif method == 'MultiMAP':
            result = anndata.read_h5ad(os.path.join(result_dir, "MultiMAP.h5ad"))
            adata.obsm[method] = result.X

        elif method == 'Seurat':
            result = anndata.read_h5ad(os.path.join(result_dir, "Seurat.h5ad"))
            adata.obsm[method] = np.array(result.X.todense())

        elif method == 'uniPort':
                result = anndata.read_h5ad(os.path.join(result_dir, "uniPort.h5ad"))
                emb = np.array(result.X.todense())
                adata.obsm[method] = np.vstack((emb[num_other:, :], emb[0:num_other, :]))

        else:
            result = np.load(os.path.join(result_dir, method + ".npy"), allow_pickle=True).item()
            adata.obsm[method] = np.vstack((result['inte'][0], result['inte'][1]))


    return adata, anno_rna, anno_other, clu_rna, clu_other

