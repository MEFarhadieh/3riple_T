"""
Gene Expression Report — Human Heart — Jupyter Compatible
Reads from adata.X directly (32732 genes, log-normalized, no adata.raw needed)
UMAP: Plotly with Spectral_r colormap + hover tooltip
"""

import json
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from pathlib import Path

# ══════════════════════════════════════════════════════
#  ★  Edit only these lines  ★
# ══════════════════════════════════════════════════════
INPUT_PATH  = "Global_lognormalised.h5ad"
OUTPUT_PATH = "col_gene_report_heart.html"
# ══════════════════════════════════════════════════════

GENES = ["GAPDH", "ACTB", "MYH6", "DCN"]

META_COLS = [
    "cell_type", "cell_state", "region",
    "donor_type", "donor", "gender",
    "cell_or_nuclei", "modality", "facility",
]

CITATION = (
    "Kanemaru, K., et al. Spatially resolved multiomics of human cardiac niches. "
    "<em>Nature</em> <strong>619</strong>, 801&ndash;810 (2023). "
    "<a href='https://doi.org/10.1038/s41586-023-06311-1' target='_blank'>"
    "https://doi.org/10.1038/s41586-023-06311-1</a>"
)

# Spectral_r: red(high) → orange → yellow → green → blue → purple(low)
SPECTRAL_R = [
    [0.0,   "rgb(94,79,162)"],
    [0.1,   "rgb(50,136,189)"],
    [0.2,   "rgb(102,194,165)"],
    [0.3,   "rgb(171,221,164)"],
    [0.45,  "rgb(230,245,152)"],
    [0.6,   "rgb(254,224,139)"],
    [0.75,  "rgb(253,174,97)"],
    [0.88,  "rgb(244,109,67)"],
    [1.0,   "rgb(158,1,66)"],
]


def safe_dense(X):
    if scipy.sparse.issparse(X):
        return np.asarray(X.todense())
    return np.asarray(X)


def get_expr(adata, gene):
    """Read from adata.X (log-normalized, 32732 genes — heart dataset has no adata.raw)."""
    if gene in adata.var_names:
        idx = adata.var_names.get_loc(gene)
        return safe_dense(adata.X[:, idx]).flatten().astype(float)
    # Fallback: search by gene_name columns in adata.var
    for col in adata.var.columns:
        if 'gene_name' in col.lower() or 'gene_id' in col.lower():
            match = adata.var[adata.var[col] == gene]
            if len(match) > 0:
                idx = adata.var.index.get_loc(match.index[0])
                print(f"  [INFO] '{gene}' found via var column '{col}'")
                return safe_dense(adata.X[:, idx]).flatten().astype(float)
    print(f"  [WARNING] '{gene}' not found in adata.var_names or var columns — skipped.")
    return None


def top_n(expr, labels):
    df = pd.DataFrame({"expr": expr, "label": labels})
    return (
        df.groupby("label")["expr"]
        .agg(
            mean_expr="mean",
            median_expr=lambda x: float(np.median(x[x > 0])) if (x > 0).any() else 0.0,
            pct_expressing=lambda x: (x > 0).mean() * 100,
            n_cells="count",
        )
        .reset_index()
        .sort_values("mean_expr", ascending=False)
    )


def compute_all_data(adata):
    avail_meta = [c for c in META_COLS if c in adata.obs.columns]
    gene_data  = {}

    for gene in GENES:
        expr = get_expr(adata, gene)
        if expr is None:
            continue
        pct = (expr > 0).mean() * 100
        print(f"  → {gene}  |  expressing: {pct:.1f}%  |  max: {expr.max():.3f}")

        per_meta = {}
        for col in avail_meta:
            labels = adata.obs[col].astype(str).values
            per_meta[col] = top_n(expr, labels).to_dict(orient="records")

        # UMAP — sample max 60k for performance
        umap_data = None
        if "X_umap" in adata.obsm:
            umap = adata.obsm["X_umap"]
            n_sample = min(60_000, adata.n_obs)
            rng = np.random.default_rng(42)
            idx = rng.choice(adata.n_obs, n_sample, replace=False)
            ct_col = "cell_type" if "cell_type" in adata.obs.columns else \
                ("Celltypes" if "Celltypes" in adata.obs.columns else None)
            ct = adata.obs[ct_col].astype(str).values[idx] if ct_col else [""] * n_sample
            expr_s = np.round(expr[idx], 4)
            umap_data = {
                "x":        np.round(umap[idx, 0], 4).tolist(),
                "y":        np.round(umap[idx, 1], 4).tolist(),
                "expr":     expr_s.tolist(),
                "celltype": ct.tolist(),
            }

        nonzero = expr[expr > 0]
        if len(nonzero) > 0:
            counts, edges = np.histogram(nonzero, bins=40)
            hist = {"counts": counts.tolist(), "edges": np.round(edges, 4).tolist()}
        else:
            hist = {"counts": [], "edges": []}

        gene_data[gene] = {
            "pct_expressing": float(pct),
            "mean_all":       float(expr.mean()),
            "median_all":     float(np.median(expr[expr > 0])) if (expr > 0).any() else 0.0,
            "max_val":        float(expr.max()),
            "n_expressing":   int((expr > 0).sum()),
            "per_meta":       per_meta,
            "umap":           umap_data,
            "hist":           hist,
        }

    found_genes = list(gene_data.keys())
    coexpr = {}
    for g1 in found_genes:
        coexpr[g1] = {}
        e1 = get_expr(adata, g1)
        for g2 in found_genes:
            e2 = get_expr(adata, g2)
            coexpr[g1][g2] = round(float(np.corrcoef(e1, e2)[0, 1]) if g1 != g2 else 1.0, 4)

    return {
        "genes":       gene_data,
        "coexpr":      coexpr,
        "meta_cols":   avail_meta,
        "found_genes": found_genes,
        "n_obs":       int(adata.n_obs),
        "n_vars":      int(adata.n_vars),
        "citation":    CITATION,
        "spectral_r":  SPECTRAL_R,
    }


HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>Collagen Gene Family Expression — Human Heart</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.2/dist/chart.umd.min.js"></script>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<link href="https://fonts.googleapis.com/css2?family=DM+Serif+Display:ital@0;1&family=DM+Sans:opsz,wght@9..40,300;9..40,400;9..40,500;9..40,600;9..40,700&family=DM+Mono:wght@400;500&display=swap" rel="stylesheet"/>
<style>
:root{
  --bg:#f5f7fa;--card:#fff;--border:#e2e8f0;--border2:#cbd5e1;
  --text:#0f172a;--muted:#475569;--muted2:#94a3b8;
  --c1:#e05a2b;--c2:#2563eb;--c3:#16a34a;--c4:#7c3aed;
  --sh0:0 1px 3px rgba(0,0,0,.07);--sh1:0 4px 18px rgba(0,0,0,.09);
}
*{box-sizing:border-box;margin:0;padding:0}
body{background:var(--bg);color:var(--text);font-family:'DM Sans',sans-serif;min-height:100vh}

/* Header */
header{
  background:linear-gradient(112deg,#eef2ff 0%,#faf5ff 50%,#ecfdf5 100%);
  border-bottom:1px solid var(--border);
  padding:2.8rem 3.5rem 2.3rem;position:relative;overflow:hidden;
}
header::before{content:'';position:absolute;top:-80px;right:-80px;width:400px;height:400px;border-radius:50%;
  background:radial-gradient(circle,rgba(124,58,237,.07) 0%,transparent 65%);pointer-events:none}
header::after{content:'';position:absolute;bottom:-70px;left:6%;width:360px;height:360px;border-radius:50%;
  background:radial-gradient(circle,rgba(37,99,235,.06) 0%,transparent 65%);pointer-events:none}
.eyebrow{font-family:'DM Mono',monospace;font-size:.67rem;letter-spacing:.22em;
  color:var(--muted2);text-transform:uppercase;margin-bottom:.85rem}
h1{font-family:'DM Serif Display',serif;font-size:2.05rem;line-height:1.15;color:var(--text);max-width:700px}
h1 em{font-style:italic;color:#7c3aed}
.hsub{color:var(--muted);font-size:.88rem;margin-top:.6rem;max-width:720px;line-height:1.65}
.cite{margin-top:1.4rem;background:#fff;border:1px solid var(--border);border-left:4px solid var(--c2);
  border-radius:8px;padding:.9rem 1.25rem;font-size:.78rem;color:var(--muted);line-height:1.7;max-width:860px}
.cite strong{color:var(--text)}.cite a{color:var(--c2);text-decoration:none}
.cite a:hover{text-decoration:underline}
.dbadge{display:inline-flex;gap:2.5rem;margin-top:1.8rem;background:#fff;
  border:1px solid var(--border);border-radius:12px;padding:.85rem 1.8rem;box-shadow:var(--sh0)}
.bval{font-family:'DM Mono',monospace;font-size:1.4rem;font-weight:500;color:var(--text)}
.bkey{font-size:.67rem;color:var(--muted2);margin-top:.2rem;letter-spacing:.04em}

/* Tabs */
nav.tabs{display:flex;gap:.6rem;padding:1.2rem 3.5rem .85rem;border-bottom:1px solid var(--border);
  background:#fff;position:sticky;top:0;z-index:100;flex-wrap:wrap;box-shadow:0 2px 10px rgba(0,0,0,.05)}
.tb{font-family:'DM Mono',monospace;font-size:.78rem;padding:.55rem 1.35rem;border-radius:8px;
  border:1.5px solid var(--border2);background:#fff;color:var(--muted);cursor:pointer;
  transition:all .18s;letter-spacing:.04em}
.tb:hover{border-color:var(--gc);color:var(--gc)}
.tb.g1{--gc:var(--c1)}.tb.g2{--gc:var(--c2)}.tb.g3{--gc:var(--c3)}.tb.g4{--gc:var(--c4)}
.tb.active{color:#fff;font-weight:500}
.tb.active.g1{background:var(--c1);border-color:var(--c1)}
.tb.active.g2{background:var(--c2);border-color:var(--c2)}
.tb.active.g3{background:var(--c3);border-color:var(--c3)}
.tb.active.g4{background:var(--c4);border-color:var(--c4)}

/* Layout */
main{padding:2rem 3.5rem;max-width:1720px;margin:0 auto}
.gp{display:none}.gp.on{display:block}

/* Stat cards */
.sgrid{display:grid;grid-template-columns:repeat(auto-fit,minmax(158px,1fr));gap:1rem;margin-bottom:2rem}
.sc{background:var(--card);border:1.5px solid var(--border);border-top:3px solid var(--gc,#e05a2b);
  border-radius:12px;padding:1.15rem 1.35rem;box-shadow:var(--sh0);transition:box-shadow .2s}
.sc:hover{box-shadow:var(--sh1)}
.sv{font-family:'DM Mono',monospace;font-size:1.55rem;font-weight:500;color:var(--gc,#e05a2b)}
.sk{font-size:.71rem;color:var(--muted);margin-top:.3rem}
.snote{font-size:.62rem;color:var(--muted2);margin-top:.15rem;font-style:italic}

/* UMAP + Hist row */
.umap-row{display:grid;grid-template-columns:540px 1fr;gap:1.5rem;margin-bottom:1.5rem;align-items:start}
@media(max-width:1150px){.umap-row{grid-template-columns:1fr}}

/* Chart card */
.cc{background:var(--card);border:1.5px solid var(--border);border-radius:16px;
  padding:1.5rem;box-shadow:var(--sh0);margin-bottom:1.5rem;transition:box-shadow .2s}
.cc:hover{box-shadow:var(--sh1)}
.ct{font-size:.7rem;font-weight:700;color:var(--muted);text-transform:uppercase;letter-spacing:.13em;
  margin-bottom:1rem;display:flex;align-items:center;gap:.6rem}
.ct::before{content:'';width:3px;height:14px;border-radius:2px;background:var(--gc,#e05a2b);flex-shrink:0}

/* Plotly UMAP container — square */
.umap-plotly{width:100%;aspect-ratio:1/1;min-height:460px}

/* Meta selector */
.mrow{display:flex;gap:.45rem;margin-bottom:1rem;flex-wrap:wrap}
.mb{font-size:.7rem;padding:.3rem .85rem;border-radius:6px;border:1.5px solid var(--border2);
  background:#fff;color:var(--muted);cursor:pointer;transition:all .15s;font-family:'DM Sans',sans-serif}
.mb:hover{border-color:var(--gc);color:var(--text)}
.mb.on{border-color:var(--gc);color:var(--gc);
  background:color-mix(in srgb,var(--gc) 8%,white);font-weight:600}

/* Bar chart */
.bwrap{position:relative;}

/* Table */
.dt{width:100%;border-collapse:collapse;font-size:.82rem}
.dt th{color:var(--text);font-weight:700;text-align:left;padding:.65rem .85rem;
  border-bottom:2px solid var(--border2);font-size:.72rem;letter-spacing:.06em;text-transform:uppercase}
.dt td{padding:.52rem .85rem;border-bottom:1px solid var(--border);vertical-align:middle}
.dt tr:last-child td{border-bottom:none}
.dt tr:hover td{background:#f8faff}
.pb{height:5px;border-radius:3px;display:inline-block;vertical-align:middle;opacity:.8}
.rk{display:inline-flex;align-items:center;justify-content:center;width:22px;height:22px;
  border-radius:50%;font-family:'DM Mono',monospace;font-size:.62rem;font-weight:500;
  background:color-mix(in srgb,var(--gc,#e05a2b) 12%,white);color:var(--gc,#e05a2b)}

/* Section label */
.sl{font-family:'DM Mono',monospace;font-size:.62rem;letter-spacing:.2em;color:var(--muted2);
  text-transform:uppercase;margin:2rem 0 1rem;display:flex;align-items:center;gap:1rem}
.sl::after{content:'';flex:1;height:1px;background:var(--border)}

/* Heatmap */
.hmc{border-radius:8px;display:flex;align-items:center;justify-content:center;
  font-family:'DM Mono',monospace;font-size:.78rem;font-weight:500;
  transition:transform .15s;min-height:65px;min-width:115px;border:1.5px solid var(--border)}
.hmc:hover{transform:scale(1.04);box-shadow:var(--sh0)}

footer{text-align:center;padding:2rem;color:var(--muted2);font-size:.72rem;
  border-top:1px solid var(--border);margin-top:3rem}
</style>
</head>
<body>

<header>
  <div class="eyebrow">Single-Cell Transcriptomics &middot; Human Lung Atlas</div>
  <h1>Collagen Gene Family Expression<br>in <em>Human Heart</em></h1>
  <div class="hsub">
    Expression patterns of GAPDH, ACTB, MYH6, and DCN across single-cell and
    single-nucleus resolved cardiac cell types, with metadata associations. Expression values from log-normalized counts (adata.X).
  </div>
  <div class="cite"><strong>Data source:</strong> __CITATION__</div>
  <div class="dbadge">
    <div><div class="bval" id="bd-obs">—</div><div class="bkey">Cells</div></div>
    <div><div class="bval" id="bd-raw">—</div><div class="bkey">Genes (raw)</div></div>
    <div><div class="bval">4</div><div class="bkey">Genes Analyzed</div></div>
  </div>
</header>

<nav class="tabs" id="tabs"></nav>
<main id="main"></main>

<footer>
  Expression from adata.raw.X (log-normalized) &middot; RNA Expression in Human Lung
</footer>

<script>
const D = __DATA_PLACEHOLDER__;
const GC  = {GAPDH:'#e05a2b',ACTB:'#2563eb',MYH6:'#16a34a',DCN:'#7c3aed'};
const GCL = ['g1','g2','g3','g4'];

const PAL = {
  GAPDH:  ['#e05a2b','#b84520','#f28660','#c93d16','#f5a880','#a33412','#f7c4a8','#7d2509'],
  ACTB: ['#2563eb','#1749c0','#6494f5','#1038a8','#93b6f8','#0c2d8f','#c0d6fb','#071d6e'],
  MYH6:  ['#16a34a','#0d7c38','#50c878','#0a5e2a','#85dda0','#074520','#aeecc2','#042e14'],
  DCN:  ['#7c3aed','#5f1fc8','#a873f5','#4a15a8','#c9a8f8','#350f80','#e0cffc','#210860'],
};

// Spectral_r colorscale for Plotly
const SPECTRAL_R = D.spectral_r;

const charts = {};
const umapDrawn = {};

const fK = n => n >= 1000 ? (n/1000).toFixed(1)+'k' : String(n);
const f  = (n, d=1) => typeof n === 'number' ? n.toFixed(d) : n;

document.getElementById('bd-obs').textContent = fK(D.n_obs);
document.getElementById('bd-raw').textContent  = fK(D.n_vars);

const tabsEl = document.getElementById('tabs');
const mainEl = document.getElementById('main');

D.found_genes.forEach((gene, gi) => {
  const gc  = GC[gene] || '#e05a2b';
  const gcl = GCL[gi]  || 'g1';

  const btn = document.createElement('button');
  btn.className = `tb ${gcl}`;
  btn.textContent = gene;
  btn.dataset.gene = gene;
  btn.onclick = () => switchTab(gene);
  tabsEl.appendChild(btn);

  const panel = document.createElement('div');
  panel.className = 'gp';
  panel.id = `p-${gene}`;
  panel.style.setProperty('--gc', gc);
  panel.innerHTML = buildPanel(gene, D.genes[gene], gc);
  mainEl.appendChild(panel);
});

const coDiv = document.createElement('div');
coDiv.innerHTML = buildCoexpr();
mainEl.appendChild(coDiv);

if (D.found_genes.length > 0) switchTab(D.found_genes[0]);

// ── Panel HTML ──────────────────────────────────────────────────────────────
function buildPanel(gene, gd, gc) {
  const mOpts = D.meta_cols.map((col, i) =>
    `<button class="mb ${i===0?'on':''}" data-col="${col}"
      onclick="switchMeta('${gene}','${col}',this)">${col}</button>`
  ).join('');

  const umapBlock = gd.umap ? `
    <div class="cc" style="--gc:${gc}">
      <div class="ct">UMAP — Expression Level</div>
      <div class="umap-plotly" id="umap-${gene}"></div>
    </div>` : '';

  return `
  <div class="sgrid" style="--gc:${gc}">
    <div class="sc" style="--gc:${gc}">
      <div class="sv">${f(gd.pct_expressing,1)}%</div>
      <div class="sk">% Expressing Cells</div>
    </div>
    <div class="sc" style="--gc:${gc}">
      <div class="sv">${fK(gd.n_expressing)}</div>
      <div class="sk">Expressing Cells</div>
    </div>
    <div class="sc" style="--gc:${gc}">
      <div class="sv">${f(gd.mean_all,4)}</div>
      <div class="sk">Mean Expression</div>
      <div class="snote">log-normalized &middot; all cells</div>
    </div>
    <div class="sc" style="--gc:${gc}">
      <div class="sv">${f(gd.median_all,4)}</div>
      <div class="sk">Median Expression</div>
      <div class="snote">expressing cells only</div>
    </div>
    <div class="sc" style="--gc:${gc}">
      <div class="sv">${f(gd.max_val,2)}</div>
      <div class="sk">Max Expression</div>
    </div>
  </div>

  <div class="umap-row">
    ${umapBlock}
    <div class="cc" style="--gc:${gc}">
      <div class="ct">Expression Distribution (Expressing Cells Only)</div>
      <div style="position:relative;height:400px"><canvas id="h-${gene}"></canvas></div>
    </div>
  </div>

  <div class="sl">Expression by Metadata</div>
  <div class="mrow" id="ms-${gene}">${mOpts}</div>
  <div class="cc" style="--gc:${gc}">
    <div class="ct">Mean Expression &amp; % Expressing — Top 20 Groups</div>
    <div class="bwrap" id="bwrap-${gene}"><canvas id="b-${gene}"></canvas></div>
  </div>

  <div class="cc" style="--gc:${gc}">
    <div class="ct">Detailed Table — All Groups (sorted by mean expression)</div>
    <div style="overflow-x:auto"><table class="dt" id="t-${gene}"></table></div>
  </div>`;
}

// ── Tab switch ──────────────────────────────────────────────────────────────
function switchTab(gene) {
  document.querySelectorAll('.tb').forEach(b => b.classList.toggle('active', b.dataset.gene === gene));
  document.querySelectorAll('.gp').forEach(p => p.classList.toggle('on', p.id === `p-${gene}`));
  if (!charts[gene]) {
    charts[gene] = {};
    requestAnimationFrame(() => requestAnimationFrame(() => initCharts(gene)));
  } else if (!umapDrawn[gene] && D.genes[gene] && D.genes[gene].umap) {
    // re-draw if container was not visible when first initialized
    const uEl = document.getElementById(`umap-${gene}`);
    if (uEl && uEl.clientWidth > 0) {
      umapDrawn[gene] = true;
      drawPlotlyUMAP(uEl, D.genes[gene].umap, gene, GC[gene] || '#e05a2b');
    }
  }
}

// ── Init all charts for a gene ──────────────────────────────────────────────
function initCharts(gene) {
  const gd  = D.genes[gene];
  const gc  = GC[gene] || '#e05a2b';
  const pal = PAL[gene] || PAL.GAPDH;

  // ── Histogram ──
  const hEl = document.getElementById(`h-${gene}`);
  if (hEl && gd.hist.counts.length > 0) {
    charts[gene].hist = new Chart(hEl, {
      type: 'bar',
      data: {
        labels: gd.hist.edges.slice(0,-1).map(v => v.toFixed(2)),
        datasets: [{
          data: gd.hist.counts,
          backgroundColor: gd.hist.counts.map((_,i) => pal[i % pal.length] + 'bb'),
          borderColor:     gd.hist.counts.map((_,i) => pal[i % pal.length]),
          borderWidth: 1, borderRadius: 3,
        }]
      },
      options: {
        responsive: true, maintainAspectRatio: false,
        plugins: { legend: { display: false } },
        scales: {
          x: {
            ticks: { color: '#1e293b', font: { size: 11, weight: '600' }, maxTicksLimit: 10 },
            grid:  { color: '#f1f5f9' },
            title: { display: true, text: 'log-normalized expression',
                     color: '#334155', font: { size: 11, weight: '600' } },
          },
          y: {
            ticks: { color: '#1e293b', font: { size: 11, weight: '600' } },
            grid:  { color: '#f1f5f9' },
            title: { display: true, text: 'Cell Count',
                     color: '#334155', font: { size: 11, weight: '600' } },
          }
        }
      }
    });
  }

  // ── Plotly UMAP ──
  const uEl = document.getElementById(`umap-${gene}`);
  if (uEl && gd.umap && !umapDrawn[gene]) {
    umapDrawn[gene] = true;
    drawPlotlyUMAP(uEl, gd.umap, gene, gc);
  }

  // ── Bar (first meta) ──
  if (D.meta_cols.length > 0) renderBar(gene, D.meta_cols[0]);
}

// ── Plotly UMAP ─────────────────────────────────────────────────────────────
function drawPlotlyUMAP(container, umap, gene, gc) {
  const { x, y, expr, celltype } = umap;
  const maxE = Math.max(...expr) || 1;

  // Split zero and expressing for better rendering
  const xi=[], yi=[], ei=[], cti=[];
  const x0=[], y0=[];

  for (let i = 0; i < x.length; i++) {
    if (expr[i] > 0) {
      xi.push(x[i]); yi.push(y[i]);
      ei.push(expr[i]); cti.push(celltype[i]);
    } else {
      x0.push(x[i]); y0.push(y[i]);
    }
  }

  const traceZero = {
    x: x0, y: y0,
    mode: 'markers',
    type: 'scattergl',
    name: 'Not expressing',
    marker: {
      color: 'rgba(203,213,225,0.35)',
      size: 2,
      line: { width: 0 },
    },
    hoverinfo: 'skip',
    showlegend: false,
  };

  const traceExpr = {
    x: xi, y: yi,
    mode: 'markers',
    type: 'scattergl',
    name: gene,
    text: cti,
    customdata: ei,
    hovertemplate:
      '<b>%{text}</b><br>' +
      'Expression: <b>%{customdata:.4f}</b><br>' +
      'UMAP1: %{x:.3f} | UMAP2: %{y:.3f}' +
      '<extra></extra>',
    marker: {
      color: ei,
      colorscale: SPECTRAL_R,
      cmin: 0,
      cmax: maxE,
      size: 3,
      opacity: 0.85,
      line: { width: 0 },
      colorbar: {
        title: { text: 'Expression<br>(log-norm)', font: { size: 11 }, side: 'right' },
        thickness: 14,
        len: 0.75,
        tickfont: { size: 10, family: 'DM Mono, monospace' },
        outlinewidth: 0,
        bgcolor: 'rgba(255,255,255,0)',
      }
    },
    showlegend: false,
  };

  const layout = {
    margin: { l: 52, r: 20, t: 20, b: 52 },
    paper_bgcolor: '#ffffff',
    plot_bgcolor:  '#ffffff',
    width:  container.clientWidth  || 520,
    height: container.clientWidth  || 520,
    xaxis: {
      title: { text: 'UMAP1', font: { size: 12, color: '#1e293b', family: 'DM Sans' } },
      showgrid: false, zeroline: false,
      tickfont: { size: 10, color: '#334155', family: 'DM Mono, monospace' },
      linecolor: '#cbd5e1', linewidth: 1.5, mirror: true,
      scaleanchor: 'y', scaleratio: 1,
    },
    yaxis: {
      title: { text: 'UMAP2', font: { size: 12, color: '#1e293b', family: 'DM Sans' } },
      showgrid: false, zeroline: false,
      tickfont: { size: 10, color: '#334155', family: 'DM Mono, monospace' },
      linecolor: '#cbd5e1', linewidth: 1.5, mirror: true,
    },
    hoverlabel: {
      bgcolor: '#fff',
      bordercolor: '#cbd5e1',
      font: { family: 'DM Sans, sans-serif', size: 12, color: '#0f172a' },
    },
    dragmode: 'pan',
  };

  const config = {
    responsive: true,
    displayModeBar: true,
    modeBarButtonsToRemove: ['select2d','lasso2d','autoScale2d'],
    scrollZoom: true,
    toImageButtonOptions: { format: 'png', filename: `UMAP_${gene}`, scale: 2 },
  };

  Plotly.newPlot(container, [traceZero, traceExpr], layout, config);
}

// ── Bar chart ───────────────────────────────────────────────────────────────
function renderBar(gene, col) {
  const gd      = D.genes[gene];
  const pal     = PAL[gene] || PAL.GAPDH;
  const allRows = gd.per_meta[col] || [];
  const rows    = allRows.slice(0, 20);   // bar: top 20 only
  if (charts[gene].bar) charts[gene].bar.destroy();

  const wrap = document.getElementById(`bwrap-${gene}`);
  if (wrap) wrap.style.height = '480px';

  charts[gene].bar = new Chart(document.getElementById(`b-${gene}`), {
    type: 'bar',
    data: {
      labels: rows.map(r => r.label),
      datasets: [
        {
          label: 'Mean Expression',
          data: rows.map(r => +r.mean_expr.toFixed(4)),
          backgroundColor: rows.map((_, i) => pal[i % pal.length] + 'bb'),
          borderColor:     rows.map((_, i) => pal[i % pal.length]),
          borderWidth: 1.5, borderRadius: 5, yAxisID: 'y',
        },
        {
          label: '% Expressing',
          data: rows.map(r => +r.pct_expressing.toFixed(1)),
          type: 'line',
          borderColor: '#475569',
          backgroundColor: 'rgba(71,85,105,.07)',
          borderWidth: 2, pointRadius: 4,
          pointBackgroundColor: '#fff',
          pointBorderColor: '#475569',
          fill: true, yAxisID: 'y2', tension: .35,
        }
      ]
    },
    options: {
      responsive: true, maintainAspectRatio: false,
      interaction: { mode: 'index', intersect: false },
      plugins: {
        legend: { labels: { color: '#334155', font: { size: 11, weight: '500' }, boxWidth: 12 } }
      },
      scales: {
        x: {
          ticks: { color: '#1e293b', font: { size: 10, weight: '600' }, maxRotation: 45 },
          grid: { color: '#f1f5f9' },
        },
        y: {
          position: 'left',
          ticks: { color: '#1e293b', font: { size: 11, weight: '600' } },
          grid: { color: '#f1f5f9' },
          title: { display: true, text: 'Mean Expression (log-norm)',
                   color: '#334155', font: { size: 11, weight: '600' } },
        },
        y2: {
          position: 'right',
          ticks: { color: '#334155', font: { size: 11, weight: '600' } },
          grid: { display: false },
          title: { display: true, text: '% Expressing',
                   color: '#334155', font: { size: 11, weight: '600' } },
        }
      }
    }
  });
  renderTable(gene, col);
}

// ── Table ────────────────────────────────────────────────────────────────────
function renderTable(gene, col) {
  const pal     = PAL[gene] || PAL.GAPDH;
  const rows    = (D.genes[gene].per_meta[col] || []);   // ALL groups
  const maxP    = Math.max(...rows.map(r => r.pct_expressing));
  document.getElementById(`t-${gene}`).innerHTML = `
  <thead><tr>
    <th>#</th><th>Group</th>
    <th>Mean Expr.</th><th>Median Expr.</th>
    <th>% Expressing</th><th>N Cells</th>
  </tr></thead>
  <tbody>${rows.map((r, i) => `
    <tr>
      <td><span class="rk" style="--gc:${pal[i%pal.length]}">${i+1}</span></td>
      <td style="font-weight:500">${r.label}</td>
      <td style="color:${pal[i%pal.length]};font-family:'DM Mono',monospace;font-weight:600">${(+r.mean_expr).toFixed(4)}</td>
      <td style="font-family:'DM Mono',monospace;color:#475569">${(+r.median_expr).toFixed(4)}</td>
      <td>
        <span class="pb" style="width:${r.pct_expressing/maxP*90}px;background:${pal[i%pal.length]}"></span>
        <span style="font-family:'DM Mono',monospace;font-size:.75rem;margin-left:.4rem;font-weight:600">${(+r.pct_expressing).toFixed(1)}%</span>
      </td>
      <td style="color:#475569;font-family:'DM Mono',monospace;font-size:.78rem">${fK(r.n_cells)}</td>
    </tr>`).join('')}
  </tbody>`;
}

// ── Meta switch ──────────────────────────────────────────────────────────────
function switchMeta(gene, col, btn) {
  document.getElementById(`ms-${gene}`)
    .querySelectorAll('.mb').forEach(b => b.classList.toggle('on', b === btn));
  renderBar(gene, col);
}

// ── Co-expression heatmap ────────────────────────────────────────────────────
function buildCoexpr() {
  const genes = D.found_genes, co = D.coexpr;
  function cs(v, gi) {
    const gc = Object.values(GC)[gi] || '#e05a2b';
    if (v >= 1) return { bg: `color-mix(in srgb,${gc} 15%,white)`, fg: gc, brd: gc };
    const abs = Math.abs(v);
    if (v > 0) {
      const a = Math.round(abs * 20);
      return { bg: `color-mix(in srgb,${gc} ${a}%,white)`, fg: abs > .3 ? '#1e293b' : '#64748b', brd: 'var(--border)' };
    }
    return { bg: `color-mix(in srgb,#ef4444 ${Math.round(abs*15)}%,white)`, fg: abs > .3 ? '#1e293b' : '#64748b', brd: 'var(--border)' };
  }
  const hdr = genes.map(g =>
    `<th style="color:${GC[g]||'#333'};font-family:'DM Mono',monospace;font-size:.75rem;padding:.5rem .8rem;font-weight:700">${g}</th>`
  ).join('');
  const rows = genes.map((g1, gi) => {
    const cells = genes.map((g2, gj) => {
      const v = co[g1] && co[g1][g2] !== undefined ? co[g1][g2] : 0;
      const { bg, fg, brd } = cs(v, gi === gj ? gi : Math.min(gi, gj));
      return `<td><div class="hmc" style="background:${bg};color:${fg};border-color:${brd}"
        title="Pearson r = ${v}">${v.toFixed(3)}</div></td>`;
    }).join('');
    return `<tr>
      <td style="color:${GC[g1]||'#333'};font-family:'DM Mono',monospace;
        font-size:.78rem;padding:.5rem .8rem;font-weight:700;white-space:nowrap">${g1}</td>
      ${cells}</tr>`;
  }).join('');

  return `
  <div class="sl" style="margin-top:3rem">Co-Expression Analysis</div>
  <div class="cc">
    <div class="ct">Pearson Correlation Matrix — COL Genes</div>
    <p style="color:var(--muted);font-size:.8rem;margin-bottom:1.2rem;line-height:1.65">
      Computed across all cells using log-normalized expression (adata.raw.X).
      Values near 1 indicate strong co-expression; values near 0 indicate independent patterns.
    </p>
    <div style="overflow-x:auto">
      <table style="border-collapse:separate;border-spacing:8px">
        <thead><tr><th></th>${hdr}</tr></thead>
        <tbody>${rows}</tbody>
      </table>
    </div>
  </div>`;
}
</script>
</body>
</html>"""


def build_html(data: dict) -> str:
    json_str = json.dumps(data, ensure_ascii=False, allow_nan=False)
    html = HTML.replace("__DATA_PLACEHOLDER__", json_str)
    html = html.replace("__CITATION__", data["citation"])
    return html


# ══════════════════════════════════════════════════════
print(f"[1/4] Loading: {INPUT_PATH}")
adata = sc.read_h5ad(INPUT_PATH)
print(f"  → {adata.n_obs:,} cells | adata.X genes: {adata.n_vars:,}  ← using adata.X directly")

print("[2/4] Computing statistics from adata.X ...")
data = compute_all_data(adata)

print("[3/4] Building HTML ...")
html = build_html(data)

out = Path(OUTPUT_PATH)
out.write_text(html, encoding="utf-8")
print(f"[4/4] ✅ Saved: {out.resolve()}  ({out.stat().st_size/1024:.0f} KB)")
