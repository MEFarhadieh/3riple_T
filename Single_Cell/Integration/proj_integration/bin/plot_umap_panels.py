#!/usr/bin/env python
import scanpy as sc, argparse, matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--infile", required=True)
parser.add_argument("--key", default="sex")   
parser.add_argument("--outfile", required=True)
parser.add_argument("--basis", default="umap")
args = parser.parse_args()

ad = sc.read_h5ad(args.infile)
sc.pl.umap(ad, color=[args.key,"batch"], wspace=0.4, show=False)
plt.savefig(args.outfile, bbox_inches="tight")
