
library(rliger)
library(Seurat)
# library(SeuratData)
library(magrittr)
library(SeuratDisk)
# library(SeuratWrappers)


data_id = "BMCITE"
dataset_dir = "../../bmcite/GASM/" 
data_batch = 's1d2_s3d7'
fea_type = 'same_name' #'combine'
rna = LoadH5Seurat(paste0(dataset_dir, "bm.rna.", data_batch, ".h5Seurat"))
adt = LoadH5Seurat(paste0(dataset_dir, "bm.adt.", data_batch, ".h5Seurat"))

cite_graph = read.csv(paste0(dataset_dir, "graph_", fea_type, '.csv'))

adt_filter <- adt@assays[["RNA"]]@data #[cite_graph[,3],]
rna_filter <- rna@assays[["RNA"]]@data #[cite_graph[,2],]


adt_name = paste0(colnames(adt_filter), '-ADT')
rna_name = paste0(colnames(rna_filter), '-RNA')
colnames(adt_filter) = adt_name
colnames(rna_filter) = rna_name

share.data <- list(adt = adt_filter, rna = rna_filter)



int.share <- createLiger(share.data)
int.share <- normalize(int.share)
int.share <- selectGenes(int.share, unshared = TRUE, 
                         unshared.datasets = list(2), 
                         unshared.thresh= 0.4)
int.share <- scaleNotCenter(int.share)
int.share <- optimizeALS(int.share, k=10, 
                         use.unshared = TRUE)
int.share <- quantile_norm(int.share, ref_dataset= "rna")
# int.share <- louvainCluster(int.share, resolution = 0.2)
# int.share <- runUMAP(int.share, distance = 'cosine', n_neighbors = 10)

combed = t(int.share@H.norm)
rownames(combed) = paste0("dim-", 1:dim(combed)[1])
combed = CreateSeuratObject(counts = combed)

celltype <- c(adt@meta.data[["cell_type2"]], rna@meta.data[["cell_type2"]])
type <- c(rep("ADT",ncol(adt_filter)), rep("RNA",ncol(rna_filter)))

combed@meta.data$cell_type=celltype
combed@meta.data$batch=type
folder_path <- paste0("../results/LIGER/", data_id)
if (!file.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
}

SaveH5Seurat(combed, filename = paste0(folder_path, '/', data_batch, ".h5Seurat"),overwrite=T)
Convert(paste0(folder_path, '/', data_batch, ".h5Seurat"), dest = "h5ad")

