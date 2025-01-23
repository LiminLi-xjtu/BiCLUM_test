rm(list = ls())
gc()

library(Matrix)
library(BiocGenerics)
library(GenomicRanges)
library(IRanges)
library(Seurat)
library(SeuratDisk)
library(zellkonverter)
library(Signac)
library(EnsDb.Hsapiens.v86)

#' Extend
#'
#' Resize GenomicRanges upstream and or downstream.
#' From \url{https://support.bioconductor.org/p/78652/}
#'
#' @param x A range
#' @param upstream Length to extend upstream
#' @param downstream Length to extend downstream
#' @param from.midpoint Count bases from region midpoint,
#' rather than the 5' or 3' end for upstream and downstream
#' respectively.
#'
#' @importFrom GenomicRanges trim
#' @importFrom BiocGenerics start strand end width
#' @importMethodsFrom GenomicRanges strand start end width
#' @importFrom IRanges ranges IRanges "ranges<-"
#' @export
#' @concept utilities
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object
#' @examples
#' Extend(x = blacklist_hg19, upstream = 100, downstream = 100)
Extend <- function(
    x,
    upstream = 0,
    downstream = 0,
    from.midpoint = FALSE
) {
  if (any(strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- strand(x = x) == "+" | strand(x = x) == "*"
  if (from.midpoint) {
    midpoints <- start(x = x) + (width(x = x) / 2)
    new_start <- midpoints - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- midpoints + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  } else {
    new_start <- start(x = x) - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- end(x = x) + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  }
  IRanges::ranges(x = x) <- IRanges::IRanges(start = new_start, end = new_end)
  x <- trim(x = x)
  return(x)
}

find_geneact <- function(peak.df, peak = peak, annotation.file, seq.levels, upstream = 2000, downstream = 0, verbose = FALSE){
  # peak.df is the regions
  # peak = peak.df
  # # reformualte peak.df of the form "chromosome", "start", "end"
  # peak.df <- do.call(what = rbind, args = strsplit(x = peak.df, split = "_"))
  # peak.df <- as.data.frame(x = peak.df)
  # colnames(x = peak.df) <- c("chromosome", "start", "end")
  
  # peak.df -> peaks.gr
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)
  BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1
  
  # gtf stores the annotation (reference genome)
  gtf <- rtracklayer::import(con = annotation.file)
  gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = "coarse")
  if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
    GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
  }
  # gtf.genes stores the genes 
  gtf.genes <- gtf[gtf$type == "gene"]
  
  # update the regions correspond to each gtf.genes, gtf.body_prom
  gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
  
  # assign peaks.gr to nearest gene region
  gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
  # only leave the ones(regions) overlap with the gene regions(distance = 0)
  keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 
                                    0]
  peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
  gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]
  gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]
  peak.ids$gene.name <- gene.ids$gene_name
  peak.ids <- as.data.frame(x = peak.ids)
  peak.ids$peak <- peak[S4Vectors::queryHits(x = keep.overlaps)]
  # new annotations should include peaks and corresponding gene.name
  annotations <- peak.ids[, c("peak", "gene.name")]
  
  return(annotations)
}



# upstream region size (base-pair)
upstream <- 2000
# downstream region size (base-pair)
downstream <- 0



#############################################################################
# hyper-parameters
species <- "Human"
data_id = "BMMC_s1d1"

# path <- "../../PBMC2/"

path <- "../../BMMC/s1d1/"


sce <- readH5AD(paste0(path,"GASM/rna.h5ad"))
data.rna = sce@assays@data@listData[["X"]]

sce2 <- readH5AD(paste0(path,"GASM/atac.h5ad"))
data.atac = sce2@assays@data@listData[["X"]]
peak_regions = rownames(data.atac)


rna <- CreateSeuratObject(counts = data.rna, project = "rna")
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna,nfeatures = 1000)
ntop_gene = rna@assays[["RNA"]]@var.features



library(stringr)
str1 = str_split(peak_regions, "-", simplify = T)
# str2 = str_split(str1[,2], "-", simplify = T)
regions = data.frame("chromosome"=str1[,1], "start"=str1[,2], "end"=str1[,3])

peak_regions = paste(str1[,1],paste(str1[,2], str1[,3], sep = '-'), sep = ':')

#############################################################################

if(species == "Mouse"){
  A = find_geneact(peak.df = regions, peak = peak_regions, annotation.file = "Mus_musculus.GRCm39.109.gtf", 
                   seq.levels = c(1:19, "X", "Y"), upstream = upstream, downstream = downstream, verbose = TRUE)
} else if(species == "Human"){
  A = find_geneact(peak.df = regions, peak = peak_regions, annotation.file = "Homo_sapiens.GRCh38.109.gtf", 
                   seq.levels = c(1:22, "X", "Y"), upstream = upstream, downstream = downstream, verbose = TRUE)
} else{
  stop("species can only be Human or Mouse")
}

A_old = A
peaks = A$peak

tmp = match(unique(A$gene.name), ntop_gene)
tmp = na.omit(tmp)
gene_names = ntop_gene[tmp]
idx = A$gene.name %in% gene_names
A = A[idx,]
peaks = peaks[idx]
peak_names = unique(peaks)

i = match(peaks, peak_names)
j = match(A$gene.name, gene_names)

str1 = str_split(peak_names, ":", simplify = T)
peak_names = paste(str1[,1],str1[,2],sep = '-')
reg = sparseMatrix(i, j, x=1)
rownames(reg) = peak_names 
colnames(reg) = gene_names

scgeneactivity = CreateSeuratObject(counts = reg)
SaveH5Seurat(scgeneactivity, filename = paste0('../data/', data_id, '.h5Seurat'), overwrite=T)
Convert(paste0('../data/', data_id, '.h5Seurat'), dest = "h5ad")
