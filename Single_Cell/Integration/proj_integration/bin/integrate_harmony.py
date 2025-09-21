#!/usr/bin/env python
import scanpy as sc, argparse
import harmonypy as hm
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--infile",  required=True)
parser.add_argument("--outfile", required=True)
parser.add_argument("--batch_key", default="batch")
parser.add_argument("--n_hvgs", type=int, default=3000)
parser.add_argument("--n_pcs",  type=int, default=50)
parser.add_argument("--neighbors", type=int, default=30)
args = parser.parse_args()

ad = sc.read_h5ad(args.infile)
sc.pp.highly_variable_genes(ad, n_top_genes=args.n_hvgs, flavor="seurat_v3", subset=True)
sc.pp.scale(ad, max_value=10)
sc.tl.pca(ad, n_comps=args.n_pcs)
ho = hm.run_harmony(ad.obsm["X_pca"], ad.obs, args.batch_key)
ad.obsm["X_harmony"] = ho.Z_corr.T
sc.pp.neighbors(ad, use_rep="X_harmony", n_neighbors=args.neighbors)
sc.tl.umap(ad)
sc.tl.leiden(ad, key_added="leiden_harmony")
ad.write(args.outfile)
