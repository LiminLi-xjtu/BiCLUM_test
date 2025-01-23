import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
import numpy as np
import os

# scglue.plot.set_publication_params()
# rcParams["figure.figsize"] = (4, 4)


# load data


def processing(rna, atac, data_id, species):
    rna.X, rna.X.data
    rna.X = rna.X.astype(int)
    rna.layers["counts"] = rna.X.copy()

    sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")

    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    sc.pp.scale(rna)
    sc.tl.pca(rna, n_comps=100, svd_solver="auto")
    # sc.pp.neighbors(rna, metric="cosine")
    # sc.tl.umap(rna)
    # sc.pl.umap(rna, color="cell_type")

    atac.X, atac.X.data
    scglue.data.lsi(atac, n_components=100, n_iter=15)
    # sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
    # sc.tl.umap(atac)
    # sc.pl.umap(atac, color="cell_type")
    if species == "Human":
        gtf = "gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz"
    else:
        gtf = "gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz"

    scglue.data.get_gene_annotation(
        rna,
        gtf=gtf,
        gtf_by="gene_name"
    )

    rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()
    tmp = rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]]
    idx = np.array(~tmp.isnull().any(axis=1))
    rna = rna[:, idx]




    atac.var_names[:5]
    split = atac.var_names.str.split(r"[:-]")
    atac.var["chrom"] = split.map(lambda x: x[0])
    atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
    atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
    atac.var.head()

    guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
    guidance

    scglue.graph.check_graph(guidance, [rna, atac])

    atac.var.head()

    path = './processed_data/'
    if not os.path.exists(path):
        os.makedirs(path)

    atac.raw.var.rename(columns={'_index': 'index'}, inplace=True)  # 将索引列名更改为 'index'
    atac.write(os.path.join(path, data_id + "-atac-pp.h5ad"), compression="gzip")

    rna.var['artif_dupl'] = rna.var['artif_dupl'].astype(str)
    rna.raw.var.rename(columns={'_index': 'index'}, inplace=True)
    rna.write(os.path.join(path, data_id + "-rna-pp.h5ad"), compression="gzip")

    nx.write_graphml(guidance, os.path.join(path, data_id + "-guidance.graphml.gz"))


