import scvi
import scanpy as sc
from scvi.model import MULTIVI
import scvi
import os



data_id = 'PBMC_paired'
species = 'Human' 
adata_rna = sc.read_h5ad('../data/'+ data_id +'/GASM/rna.h5ad')
adata_atac = sc.read_h5ad('../data/'+ data_id +'/GASM/atac.h5ad')

adata_rna.layers["rna"] = adata_rna.X.copy()
# adata_rna.layers["atac"] = None
adata_rna.obs["modality"] = "rna"

adata_atac.layers["atac"] = adata_atac.X.copy()
# adata_atac.layers["rna"] = None
adata_atac.obs["modality"] = "atac"



adata = adata_rna.concatenate(adata_atac, batch_key="batch", join="outer")
print(adata.obs["batch"].value_counts())


MULTIVI.setup_anndata(
    adata,
    layer=None,   
    batch_key="batch"
)

# scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")


model = MULTIVI(
    adata,
    n_genes=adata_rna.shape[1],
    n_regions=adata_atac.shape[1],
)


model.train()


latent = model.get_latent_representation()
MULTIVI_LATENT_KEY = "X_multivi"
adata.obsm[MULTIVI_LATENT_KEY] = latent

rna_res = adata[:adata_rna.shape[0], :]
atac_res = adata[adata_rna.shape[0]:, :]

inte1 = sc.AnnData(rna_res.obsm['X_multivi'])
inte2 = sc.AnnData(atac_res.obsm['X_multivi'])
inte1.obs_names = rna_res.obs_names
inte2.obs_names = atac_res.obs_names
inte1.obs['cell_type'] = rna_res.obs['cell_type']
inte2.obs['cell_type'] = atac_res.obs['cell_type']
inte1.obs['batch'] = rna_res.obs['batch']
inte2.obs['batch'] = atac_res.obs['batch']


path = './results/'+ data_id
if not os.path.exists(path):
    os.makedirs(path)

inte1.write(os.path.join(path, "multivi_rna.h5ad"), compression="gzip")
inte2.write(os.path.join(path, "multivi_atac.h5ad"), compression="gzip")




