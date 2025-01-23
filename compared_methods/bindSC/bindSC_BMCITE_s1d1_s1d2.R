# install_github('KChen-lab/bindSC')
# rmarkdown::render('bindsc_linux.Rmd',output_file='ArchR.html')
library(zellkonverter)
library(bindSC)
library(Seurat)
library(Signac)
library(SeuratDisk)
# library(SeuratData)
library(peakRAM)

data_id = "BMCITE"
dataset_dir = "../../data/bmcite/GASM/"
data_batch = 's1d1_s1d2'

rna = LoadH5Seurat(paste0(dataset_dir, "bm.rna.", data_batch, ".h5Seurat"))
adt = LoadH5Seurat(paste0(dataset_dir, "bm.adt.", data_batch, ".h5Seurat"))

cite_graph = read.csv(paste0(dataset_dir, 'graph.csv'))
var_rna = rna@assays[["RNA"]]@var.features

y <- rna@assays[["RNA"]]@data[var_rna, ]
# y <- bm@assays[["RNA"]]@data[bm@assays[["RNA"]]@var.features, ]
x <- adt@assays[["RNA"]]@data[cite_graph[,3],]
z0 <- rna@assays[["RNA"]]@data[cite_graph[,2],]

dim(x)
dim(y)
dim(z0)
x.clst <- rna$cell_type2
y.clst <- adt$cell_type2


X <- x
Z0 <-  z0

celltype <- c(as.character(x.clst), as.character(y.clst))
type <- c(rep("ADT",ncol(X)), rep("RNA",ncol(y)))


run = FALSE
if(run){
  down_sample <- seq(1,N,10)
  paraSel <- BiCCA_para_opt( X = x[, down_sample] ,
                             Y = y$dt1[, down_sample], 
                             Z0 =z0[, down_sample], 
                             tp="out",
                             X.clst = x.clst[down_sample],
                             Y.clst = y.clst[down_sample],
                             alpha.lst = seq(0,1,0.1), 
                             K.lst = c(10),
                             lambda.lst = seq(0,1,0.1),
                             num.iteration = 50,
                             tolerance = 0.01,
                             save = TRUE,
                             block.size = 1000
  )
  p1 <- paraSel_plot(paraSel)
  p1
}



res <- BiCCA( X = x ,
              Y =  y, 
              Z0 =z0, 
              X.clst = x.clst,
              Y.clst = y.clst,
              alpha = 0.1, 
              lambda = 0.7,
              K = 15,
              temp.path  = "out",
              num.iteration = 50,
              tolerance = 0.01,
              save = TRUE,
              parameter.optimize = FALSE,
              block.size = 0)

rownames(res$r) = paste0(rownames(res$r), "-2")
cobed = rbind(res$u, res$r)
cobed = t(cobed)
rownames(cobed) = paste0("dim-", 1:dim(cobed)[1])
result = list(cobed, celltype, type)
saveRDS(result, paste0('bindSC_', data_id, '_', data_batch, '_', fea_type, '.rds'))

# cobed = readRDS('results/bindsc/bmcite/bindSC_bmcite.rds')
cobed = CreateSeuratObject(counts = cobed)
cobed@meta.data[["cell_type"]]=c(x.clst, y.clst)
cobed@meta.data[["batch"]]=type
folder_path <- paste0("../results/bindsc/", data_id)
if (!file.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
}
SaveH5Seurat(cobed, filename = paste0(folder_path, '/', data_batch, '_', fea_type, ".h5Seurat"),overwrite=T)
Convert(paste0(folder_path, '/', data_batch, '_', fea_type, ".h5Seurat"), dest = "h5ad")


