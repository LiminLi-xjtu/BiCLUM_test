
from pandas import value_counts
import scanpy as sc
import scipy
import numpy as np 
import warnings
warnings.filterwarnings('ignore')
from sklearn.preprocessing import normalize
import scipy.sparse as sps
import sklearn

def tfidf_data(X):
    r"""
    TF-IDF normalization (following the Seurat v3 approach)
        TF-IDF normalized matrix
    """
    idf = X.shape[0] / X.sum(axis=0)
    if scipy.sparse.issparse(X):
        tf = X.multiply(1 / X.sum(axis=1))
        return tf.multiply(idf)
    else:
        tf = X / X.sum(axis=1, keepdims=True)
        return tf * idf

def atac_lsi(adata, n_components = 20,use_highly_variable = None, random_state=0):
    # Keep deterministic as the default behavior
    if use_highly_variable is None:
        use_highly_variable = "highly_variable" in adata.var
    adata_use = adata[:, adata.var["highly_variable"]] if use_highly_variable else adata
    X = tfidf_data(adata_use.X)
    X_norm = normalize(X, norm="l1")
    X_norm = np.log1p(X_norm * 1e4)
    X_lsi = sklearn.utils.extmath.randomized_svd(X_norm, n_components)[0]
    X_lsi -= X_lsi.mean(axis=1, keepdims=True)
    X_lsi /= X_lsi.std(axis=1, ddof=1, keepdims=True)
    adata.obsm["X_lsi"] = X_lsi

    
def rna_pca(adata, n_components = 20):
    # sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=True)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, zero_center=True, max_value=10)
    sc.tl.pca(adata, n_comps=n_components, zero_center=True, svd_solver='arpack')

