#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(SeuratDisk); library(optparse) })

if (!requireNamespace("rliger", quietly = TRUE)) {
  install.packages("remotes")
  remotes::install_github("MacoskoLab/liger")
}
suppressPackageStartupMessages({ library(rliger); library(Seurat); library(dplyr) })

option_list <- list(
  make_option("--infile",  type="character"),
  make_option("--outfile", type="character"),
  make_option("--batch_key", type="character", default="batch"),
  make_option("--neighbors", type="integer", default=30)
)
opt <- parse_args(OptionParser(option_list=option_list))

obj <- SeuratDisk::LoadH5AD(opt$infile)
# split by batch → list of matrices
lst <- lapply(split(colnames(obj), obj@meta.data[[opt$batch_key]]), function(cols) obj@assays$RNA@data[, cols])
lig <- createLiger(lst)
lig <- normalize(lig); lig <- selectGenes(lig); lig <- scaleNotCenter(lig)
lig <- optimizeALS(lig, k=30); lig <- quantile_norm(lig)
emb <- lig@H.norm[[1]] 
# for UMAP with Seurat compatibility
obj[["liger"]] <- CreateDimReducObject(emb, key="liger_", assay=DefaultAssay(obj))
obj <- FindNeighbors(obj, reduction="liger") %>% RunUMAP(reduction="liger") %>% FindClusters(resolution=0.8)
SeuratDisk::SaveH5Seurat(obj, filename="tmp.h5seurat", overwrite=TRUE)
SeuratDisk::Convert("tmp.h5seurat", dest="h5ad", overwrite=TRUE)
file.rename("tmp.h5ad", opt$outfile)
