import anndata
import scanpy as sc
import sccross
import os


# dataset_name = 'PBMC_unpaired'
# species = 'Human' 
# rna = anndata.read_h5ad('../../data/PBMC_unpaired/GASM/rna.h5ad')
# atac = anndata.read_h5ad('../../data/PBMC_unpaired/GASM/atac.h5ad')


# dataset_name = 'PBMC_paired'
# species = 'Human' 
# rna = anndata.read_h5ad('../../data/PBMC_paired/GASM/rna.h5ad')
# atac = anndata.read_h5ad('../../data/PBMC_paired/GASM/atac.h5ad')


# dataset_name = 'Kidney'
# species = 'Human' 
# rna = anndata.read_h5ad('../../data/Kidney/GASM/rna.h5ad')
# atac = anndata.read_h5ad('../../data/Kidney/GASM/atac.h5ad')


# dataset_name = 'BMMC_paired'
# species = 'Human' 
# rna = anndata.read_h5ad('../../data/BMMC_s1d1/GASM/rna.h5ad')
# atac = anndata.read_h5ad('../../data/BMMC_s1d1/GASM/atac.h5ad')

dataset_name = 'BMMC_unpaired'
species = 'Human' 
rna = anndata.read_h5ad('../data/BMMC_s1d1/GASM/rna.h5ad')
atac = anndata.read_h5ad('../data/BMMC_s1d2/GASM/atac.h5ad')



rna.X = rna.layers["counts"].copy()
sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")

sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")
sc.pp.neighbors(rna, metric="cosine")
sc.tl.umap(rna)

sccross.data.lsi(atac, n_components=100, n_iter=15)
sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
sc.tl.umap(atac)

atac2rna = sccross.data.geneActivity(atac, gtf_file='./reference/gencode.v41.annotation.gtf')

path = "data/"+dataset_name+ "/"
if not os.path.exists(path):
    os.makedirs(path)


if rna.raw is not None and "_index" in rna.raw.var.columns:
        del rna.raw.var["_index"]

if atac.raw is not None and "_index" in atac.raw.var.columns:
        del atac.raw.var["_index"]

if atac2rna.raw is not None and "_index" in atac2rna.raw.var.columns:
        del atac2rna.raw.var["_index"]


atac2rna.write(path+"atac_preprocessed.h5ad", compression="gzip")

rna.write(path+"rna_preprocessed.h5ad", compression="gzip")


