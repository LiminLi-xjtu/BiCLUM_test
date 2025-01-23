
library(Seurat)
library(Signac)
library(ggplot2)
library(cowplot)
library("SeuratDisk")
# library(SeuratData)
library(magrittr)

rm(list=ls())


data_id = "BMCITE"
dataset_dir = '../../data/bmcite/GASM/'
data_batch = 's1d2_s3d7'


cite_graph = read.csv(paste0(dataset_dir, "graph.csv"))


cite.rna = LoadH5Seurat(paste0(dataset_dir, "bm.rna.", data_batch, ".h5Seurat"))
adt = LoadH5Seurat(paste0(dataset_dir, "bm.adt.", data_batch, ".h5Seurat"))

adt_data = adt@assays[["RNA"]]@counts
idx = match(cite_graph[,3], rownames(adt))
rownames(adt_data)[idx] = cite_graph[,2]
cite.adt <- CreateSeuratObject(counts = adt_data, assay = "ADT")
VariableFeatures(cite.adt) <- rownames(adt_data)
cite.adt <- NormalizeData(cite.adt, normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>% RunPCA(reduction.name = 'apca')


celltype <- c(cite.adt@meta.data[["cell_type2"]], cite.rna@meta.data[["cell_type2"]])
type <- c(rep("ADT",ncol(cite.adt)), rep("RNA",ncol(cite.rna)))


gene.use = cite_graph[,2]
# transfer.anchors <- FindTransferAnchors(reference = cite.rna, query = cite.adt, features = VariableFeatures(object = cite.rna), reference.assay = "RNA", query.assay = "ADT", reduction = "cca")
transfer.anchors <- FindTransferAnchors(reference = cite.rna,
                                        query = cite.adt,
                                        features = gene.use,
                                        dims = seq(1,15,1),
                                        reference.assay = "RNA",
                                        query.assay = "ADT",
                                        reduction = "cca")

refdata <- GetAssayData(cite.rna, assay = "RNA", slot = "data")[gene.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = cite.adt[["apca"]],
                           dims = 2:20)
cite.adt[["RNA"]] <- imputation

coembed <- merge(x = cite.rna, y = cite.adt)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = gene.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = gene.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:20)


folder_path <- paste0("results/Seurat/", data_id)
if (!file.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
}



result = list(coembed, celltype, type)


coembed = t(result[[1]]@reductions[["pca"]]@cell.embeddings)
cell_type = result[[2]]
type = result[[3]]
coembed <- CreateSeuratObject(counts = coembed, assay = "RNA")
coembed@meta.data$cell_type=cell_type
coembed@meta.data$batch=type


SaveH5Seurat(coembed, filename = paste0(folder_path, '/', data_batch, '_', fea_type, ".h5Seurat"),overwrite=T)
Convert(paste0(folder_path, '/', data_batch, '_', fea_type, ".h5Seurat"), dest = "h5ad")

##########################################################################


