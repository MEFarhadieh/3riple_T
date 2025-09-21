#!/usr/bin/env python
import argparse, scanpy as sc, anndata as ad, pandas as pd, numpy as np
from sklearn.metrics import silhouette_score, roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
import scib

parser = argparse.ArgumentParser()
parser.add_argument("--infile", required=True)
parser.add_argument("--embed_key", required=True)   
parser.add_argument("--batch_key", default="batch")
parser.add_argument("--sex_key", default="sex")
parser.add_argument("--label_key", default=None)    
parser.add_argument("--outfile", required=True)
args = parser.parse_args()

ad = sc.read_h5ad(args.infile)
emb = ad.obsm[args.embed_key]
batch = ad.obs[args.batch_key].astype(str)
sex   = ad.obs[args.sex_key].astype(str)

# --- scIB metrics ---
scores = {}
:
try:
    scores["iLISI"] = scib.me.lisi_graph(ad, batch_key=args.batch_key, type_="iLISI", use_rep=args.embed_key)["iLISI"].mean()
except Exception: scores["iLISI"]=np.nan
try:
    scores["kBET"]  = scib.me.kBET(ad, batch_key=args.batch_key, embed=args.embed_key)
except Exception: scores["kBET"]=np.nan
try:
    scores["graph_conn"] = scib.me.graph_connectivity(ad, label_key=args.label_key) if args.label_key else np.nan
except Exception: scores["graph_conn"]=np.nan

# --- Preserving gender differences (Bio signal preservation) ---
# 1) silhouette based on sex
try:
    scores["ASW_sex"] = silhouette_score(emb, sex)
except Exception: scores["ASW_sex"]=np.nan
# 2) AUC of sex prediction from embedding
try:
    le = (sex.astype("category").cat.codes).values
    aucs=[]
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)
    for tr, te in skf.split(emb, le):
        clf = LogisticRegression(max_iter=200).fit(emb[tr], le[tr])
        proba = clf.predict_proba(emb[te])[:,1]
        aucs.append(roc_auc_score(le[te], proba))
    scores["AUC_sex"] = float(np.mean(aucs))
except Exception: scores["AUC_sex"]=np.nan
# 3) Variance explained (R2) for sex (with logistic regression as an approximation)
try:
    clf = LogisticRegression(max_iter=200).fit(emb, le)
    scores["R2_sex_pseudo"] = clf.score(emb, le)
except Exception: scores["R2_sex_pseudo"]=np.nan

pd.DataFrame([scores]).to_csv(args.outfile, index=False)
