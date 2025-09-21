#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(dplyr)
  library(optparse); library(SeuratDisk); library(Matrix)
})

option_list <- list(
  make_option("--infile",  type="character"),
  make_option("--outfile", type="character"),
  make_option("--batch_key", type="character", default="batch"),
  make_option("--n_hvgs", type="integer", default=3000),
  make_option("--n_pcs",  type="integer", default=50),
  make_option("--neighbors", type="integer", default=30)
)
opt <- parse_args(OptionParser(option_list=option_list))

library(SeuratDisk)
ad <- SeuratDisk::LoadH5AD(opt$infile)    
obj <- ad
# split by batch
objs <- SplitObject(obj, split.by = opt$batch_key)
objs <- lapply(objs, function(x){ x <- NormalizeData(x); x <- FindVariableFeatures(x, nfeatures=opt$n_hvgs); x})
features <- SelectIntegrationFeatures(objs, nfeatures=opt$n_hvgs)
objs <- lapply(objs, function(x){ x <- ScaleData(x, features=features); x <- RunPCA(x, features=features); x})
anchors <- FindIntegrationAnchors(objs, anchor.features=features, reduction="rpca")
combined <- IntegrateData(anchorset=anchors)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined) %>% RunPCA(npcs=opt$n_pcs) %>% FindNeighbors(dims=1:opt$n_pcs, k.param=opt$neighbors) %>% RunUMAP(dims=1:opt$n_pcs)
combined <- FindClusters(combined, resolution=0.8)
# output h5ad
SaveH5Seurat(combined, filename="tmp.h5seurat", overwrite=TRUE)
Convert("tmp.h5seurat", dest="h5ad", overwrite=TRUE)
file.rename("tmp.h5ad", opt$outfile)
