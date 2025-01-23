from itertools import chain

import anndata as ad
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
from matplotlib import rcParams
import os

# scglue.plot.set_publication_params()
# rcParams["figure.figsize"] = (4, 4)
path = './processed_data/'

def training(data_id):
    rna = ad.read_h5ad(os.path.join(path, data_id + "-rna-pp.h5ad"))
    atac = ad.read_h5ad(os.path.join(path, data_id + "-atac-pp.h5ad"))
    guidance = nx.read_graphml(os.path.join(path, data_id + "-guidance.graphml.gz"))

    # features = rna.var.query("highly_variable").index.to_numpy().tolist()

    scglue.models.configure_dataset(
        rna, "NB", use_highly_variable=True,
        use_layer="counts", use_rep="X_pca"
    )

    scglue.models.configure_dataset(
        atac, "NB", use_highly_variable=True,
        use_rep="X_lsi"
    )

    guidance_hvf = guidance.subgraph(chain(
        rna.var.query("highly_variable").index,
        atac.var.query("highly_variable").index
    )).copy()

    glue = scglue.models.fit_SCGLUE(
        {"rna": rna, "atac": atac}, guidance_hvf,
        fit_kws={"directory": "glue"}
    )

    dx = scglue.models.integration_consistency(
        glue, {"rna": rna, "atac": atac}, guidance_hvf
    )
    dx

    # _ = sns.lineplot(x="n_meta", y="consistency", data=dx).axhline(y=0.05, c="darkred", ls="--")

    rna.obsm["X_glue"] = glue.encode_data("rna", rna)
    atac.obsm["X_glue"] = glue.encode_data("atac", atac)
    # combined = ad.concat([rna, atac])
    # sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
    # sc.tl.umap(combined)
    # sc.pl.umap(combined, color=["cell_type", "batch"], wspace=0.65)

    return rna, atac
