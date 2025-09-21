#!/usr/bin/env python
import os, sys, scanpy as sc, pandas as pd, anndata as ad
from pathlib import Path

def read_input(path):
    p = Path(path)
    if p.suffix.lower()==".h5ad":
        return sc.read_h5ad(path)
    if p.suffix.lower()==".rds":
        raise RuntimeError(f"RDS found: {path}. Please convert to h5seurat and then h5ad with SeuratDisk.")
    if p.is_dir():
        return sc.read_10x_mtx(path, var_names='gene_symbols', cache=True)
    raise RuntimeError(f"Unknown format: {path}")

meta = pd.read_csv("meta/samplesheet.csv")
adatas = []
for _, row in meta.iterrows():
    adata = read_input(os.path.expanduser(row["file"]))
    adata.obs["sample"] = row["sample"]
    adata.obs["batch"]  = row["batch"]
    adata.obs["sex"]    = row["sex"]
    if "label" in row and pd.notnull(row["label"]):
        adata.obs["label"] = str(row["label"])
    adata.var_names_make_unique()
    adatas.append(adata)

adata = ad.concat(adatas, join="outer", label="batch", index_unique="-")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write("data/combined_raw.h5ad")
print("Wrote data/combined_raw.h5ad")
