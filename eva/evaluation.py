

from metrics import *

def evaluate(adata, anno_rna, anno_other, paired):

    eva_metrics = {}
    obsm_names = adata.obsm_keys()
    for i in range(len(obsm_names)):
        method = obsm_names[i]
        print(method)
        tmp = all_metrics(adata.obsm[method], np.array(adata.obs['cluster'].array), anno_rna, anno_other,
                          np.array(adata.obs['omic_id'].array), paired=paired)
        eva_metrics[method] = tmp

    eva_metrics = pd.DataFrame(eva_metrics)
    eva_metrics.index = ['omics_mixing', 'biology_conservation', 'transfer_accuracy', 'foscttm']

    return eva_metrics

