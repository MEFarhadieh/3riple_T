import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import issparse
from scipy.stats import median_abs_deviation
import scipy.stats
from statsmodels.stats.multitest import multipletests
import scrublet as scr
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# R integration
import logging
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro

rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

# Set plotting parameters
sc.settings.set_figure_params(dpi=150, frameon=False, figsize=(8, 6))
sc.settings.verbosity = 1


def is_outlier(adata, metric: str, nmads: int):
    """Detect outliers using MAD (Median Absolute Deviation)"""
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def test_outlier(x, upper_mad_only=True):
    """Test for outliers using MAD with statistical testing"""
    med = np.median(x)
    if upper_mad_only:
        mad = np.median(x[x>med] - med) * 1.4826
    else:
        mad = np.median(np.abs(x - med)) * 1.4826
    pvals = 1 - scipy.stats.norm.cdf(x, loc=med, scale=mad)
    bh_pvals = multipletests(pvals, method='fdr_bh')[1]
    return pvals, bh_pvals


def create_qc_plots(adata, title_prefix=""):
    """Create comprehensive QC plots"""
    plots = {}
    
    # Plot 1: Distributions
    sns.set_style("ticks")
    fig1, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=False)
    sns.histplot(adata.obs['n_genes_by_counts'], ax=ax1, kde=True, bins=100)
    ax1.set_title(f'{title_prefix}Genes per cell')
    ax1.set_xlabel('Number of genes')
    
    sns.histplot(adata.obs['total_counts'], ax=ax2, kde=True, bins=100)
    ax2.set_title(f'{title_prefix}UMI counts per cell')
    ax2.set_xlabel('Total counts')
    
    sns.histplot(adata.obs['pct_counts_MT'], ax=ax3, kde=True, bins=100)
    ax3.set_title(f'{title_prefix}MT% per cell')
    ax3.set_xlabel('MT%')
    
    fig1.text(-0.01, 0.5, 'Frequency', ha='center', va='center', 
              rotation='vertical', size='x-large')
    fig1.tight_layout()
    plots['distributions'] = fig1
    
    # Plot 2: Scatter plot
    sns.set_style("ticks")
    fig2, ax = plt.subplots(figsize=(8, 6), dpi=150)
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", 
                  color="pct_counts_MT", color_map="flare", ax=ax, show=False)
    ax.set_title(f"{title_prefix}Total counts vs Genes")
    plt.tight_layout()
    plots['scatter'] = fig2
    
    # Plot 3: Highest expressed genes - capture the figure properly
    sns.set_style("whitegrid")
    fig3, ax3 = plt.subplots(figsize=(8, 6), dpi=150)
    sc.pl.highest_expr_genes(adata, n_top=20, ax=ax3, show=False)
    ax3.set_title(f"{title_prefix}Top 20 highest expressed genes")
    plt.tight_layout()
    plots['highest_expr'] = fig3
    
    return plots


def run_soupx(adata, raw_matrix_path, sample_name, output_dir):
    """Run SoupX contamination removal using R"""
    print("  Running SoupX contamination removal...")
    
    # Preprocessing for clustering
    print("    Preprocessing for clustering...")
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=min(50, min(adata_pp.n_obs, adata_pp.n_vars) - 1))
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="soupx_groups")
    
    # Get cluster assignments - convert to proper format
    soupx_groups = adata_pp.obs["soupx_groups"].astype(str).values
    n_clusters = len(np.unique(soupx_groups))
    print(f"    Found {n_clusters} clusters for SoupX")
    
    del adata_pp
    
    # Prepare data for R
    cells = np.array(adata.obs_names.tolist())
    genes = np.array(adata.var_names.tolist())
    data = adata.X.copy()
    
    if not issparse(data):
        from scipy.sparse import csr_matrix
        data = csr_matrix(data)
    
    # Load raw matrix
    print(f"    Loading raw matrix from: {raw_matrix_path}")
    if not os.path.exists(raw_matrix_path):
        raise FileNotFoundError(f"Raw matrix not found: {raw_matrix_path}")
    
    adata_raw = sc.read_10x_h5(filename=raw_matrix_path)
    adata_raw.var_names_make_unique()
    data_tod = adata_raw.X.copy()
    
    if not issparse(data_tod):
        from scipy.sparse import csr_matrix
        data_tod = csr_matrix(data_tod)
    
    print(f"    Filtered matrix: {data.shape}, Raw matrix: {data_tod.shape}")
    del adata_raw
    
    # Transfer to R - use proper conversion
    print("    Transferring data to R...")
    ro.globalenv['data'] = data.T  # Transpose for R (genes x cells)
    ro.globalenv['data_tod'] = data_tod.T
    ro.globalenv['genes'] = ro.StrVector(genes)
    ro.globalenv['cells'] = ro.StrVector(cells)
    ro.globalenv['soupx_groups'] = ro.StrVector(soupx_groups)
    
    # Run SoupX in R
    print("    Running SoupX in R...")
    ro.r(f'''
    suppressMessages(library(SoupX))
    
    # Ensure proper format
    data <- as(data, "sparseMatrix")
    data_tod <- as(data_tod, "sparseMatrix")
    
    # Set names
    rownames(data) <- genes
    colnames(data) <- cells
    
    # Create SoupChannel
    sc <- SoupChannel(data_tod, data, calcSoupProfile = FALSE)
    
    # Set soup profile
    soupProf <- data.frame(
        row.names = rownames(data), 
        est = rowSums(data)/sum(data), 
        counts = rowSums(data)
    )
    sc <- setSoupProfile(sc, soupProf)
    
    # Set clusters
    sc <- setClusters(sc, soupx_groups)
    
    # Estimate contamination and save plot
    png("{output_dir}/soupx_contamination_plot.png", width=800, height=600)
    sc <- autoEstCont(sc, doPlot=TRUE)
    dev.off()
    
    # Get contamination fraction
    rho_values <- sc$metaData$rho
    contamination_fraction <- mean(rho_values, na.rm=TRUE)
    
    # Adjust counts
    out <- adjustCounts(sc, roundToInt = TRUE)
    
    # Convert back to matrix format for transfer
    out <- as.matrix(out)
    
    cat("SoupX completed successfully\n")
    cat("Contamination fraction:", contamination_fraction, "\n")
    cat("Output dimensions:", dim(out), "\n")
    ''')
    
    # Get results back from R
    print("    Retrieving results from R...")
    try:
        out = ro.globalenv['out']
        contamination_fraction = ro.globalenv['contamination_fraction']
        
        # Convert to numpy - handle different R return types
        if hasattr(out, 'to_numpy'):
            out = out.to_numpy()
        else:
            out = np.array(out)
        
        # Extract scalar from contamination_fraction
        if hasattr(contamination_fraction, '__len__') and len(contamination_fraction) > 0:
            contamination_fraction = float(contamination_fraction[0])
        else:
            contamination_fraction = float(contamination_fraction)
        
        print(f"    Estimated contamination fraction: {contamination_fraction:.2%}")
        print(f"    Corrected counts shape: {out.shape}")
        
        # Verify dimensions
        expected_shape = (len(genes), len(cells))
        if out.shape != expected_shape:
            print(f"    Warning: Shape mismatch. Expected {expected_shape}, got {out.shape}")
            print(f"    Transposing...")
            out = out.T
        
    except Exception as e:
        print(f"    ERROR extracting SoupX results: {e}")
        print(f"    out type: {type(ro.globalenv.get('out', 'NOT FOUND'))}")
        print(f"    contamination_fraction type: {type(ro.globalenv.get('contamination_fraction', 'NOT FOUND'))}")
        raise ValueError(f"SoupX failed to return valid results: {e}")
    
    # Store in adata
    print("    Storing SoupX results...")
    adata.layers["counts"] = adata.X.copy()
    
    # Convert to sparse if needed
    if issparse(adata.X):
        from scipy.sparse import csr_matrix
        adata.layers["soupX_counts"] = csr_matrix(out.T)
    else:
        adata.layers["soupX_counts"] = out.T
    
    adata.X = adata.layers["soupX_counts"]
    adata.uns['soupx_contamination'] = contamination_fraction
    
    print("    SoupX completed successfully!")
    
    return adata


def run_scrublet(adata, resolution_function=None, output_dir=None):
    """Run Scrublet doublet detection with clustering-based refinement"""
    print(f"  Running advanced Scrublet with clustering...")
    
    old_verbosity = sc.settings.verbosity
    sc.settings.verbosity = 1
    
    if resolution_function is None:
        resolution_function = lambda x: np.maximum(np.maximum(np.log10(x)-1, 0)**2, 0.1)
    
    # Run basic Scrublet
    print("    Running Scrublet doublet detection...")
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06, random_state=42)
    
    try:
        ds, pd = scrub.scrub_doublets(verbose=False)
    except Exception as e:
        print(f"    ERROR: Scrublet failed: {e}")
        return adata, None
    
    adata.obs['scrublet_score'] = ds
    
    # Create copy for clustering
    print("    Preprocessing for clustering...")
    adata_copy = adata.copy()
    sc.pp.filter_genes(adata_copy, min_cells=3)
    sc.pp.normalize_total(adata_copy, target_sum=1e4)
    sc.pp.log1p(adata_copy)
    sc.pp.highly_variable_genes(adata_copy, flavor='seurat_v3', n_top_genes=4000)
    sc.pp.scale(adata_copy, zero_center=False)
    sc.pp.pca(adata_copy, svd_solver='arpack', zero_center=False)
    sc.pp.neighbors(adata_copy, n_pcs=30)
    sc.tl.umap(adata_copy)
    
    # Initial clustering
    print("    Performing initial clustering...")
    sc.tl.leiden(adata_copy, resolution=1)
    
    # Create UMAP plots
    print("    Creating UMAP plots...")
    fig1 = plt.figure(figsize=(16, 6), dpi=150)
    sc.pl.umap(adata_copy, color=['leiden','scrublet_score'], 
               show=False, title=['Leiden Clusters', 'Scrublet Score'])
    plt.tight_layout()
    if output_dir:
        plt.savefig(os.path.join(output_dir, 'umap_leiden_scrublet.png'), dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()
    
    # Create marker dotplots
    marker_dict = {
        "Nuclear genes": ['MALAT1', 'NEAT1', 'FTX', 'FOXP1', 'RBMS3', 'ZBTB20', 'LRMDA', 'PBX1', 'ITPR2', 'AUTS2'],
        "Mitochondrial genes": ['MT-ND1', 'MT-ND2', 'MT-CO1', 'MT-CO2', 'MT-ATP8', 'MT-ATP6', 'MT-CO3', 'MT-ND3', 'MT-ND4L', 'MT-ND4', 'MT-ND5', 'MT-ND6', 'MT-CYB'],
        "VEC": ["VWF", "ERG", "ANO2", "PTPRB", "EGFL7", "PREX2", "ADGRL4", "FLT1", "CYYR1", "GRB10", "PPP1R16B", "DOCK9", "SHANK3", "PECAM1", "PLEKHG1", "EMCN"],
        "PER": ["RGS5", "DACH1", "GUCY1A1", "ABCC9", "BGN", "NOTCH3", "PDGFRB", "FRMD3", "RNF152", "CCDC102B", "NGF"],
        "SMC": ["MYH11", "KCNAB1", "NTRK3", "CHRM3", "ACTA2", "RGS6", "DGKG", "ITGA8", "TBX2", "LMOD1", "SDK1", "GPC6", "ANTXR1", "FLNA", "CLMN", "ATP10A", "MCAM", "TAGLN", "CCDC3"],
        "AD": ["PLIN4", "PLIN1", "PDE3B", "GPAM", "PTPRS", "PPARG", "MLXIPL", "MGST1", "AQP7", "SLC19A3", "FABP4", "TPRG1", "DIRC3", "LPL", "PNPLA2", "LIPE", "ADH1B", "ADIPOQ"],
        "SC": ["XKR4", "SLC35F1", "ZNF536", "NCAM2", "GPM6B", "KIRREL3", "SORCS1", "ST6GALNAC5", "PRKCA", "GINS3", "PMP22", "ALDH1A1", "IL1RAPL2", "DOCK5", "NKAIN3"],
        "N": ["CSMD1", "SYT1", "KCNIP4", "CNTNAP2", "DLGAP1", "PTPRD", "LRRTM4", "ATRNL1", "LRP1B", "CTNND2", "KCNQ5", "NRG3", "SNTG1", "GRIA2", "RIMS2", "CSMD3"],
        "EEC": ["PCDH7", "PCDH15", "LINC02147", "LINC02388", "MYRIP", "GMDS", "ADAMTSL1", "LEPR", "CALCRL", "CGNL1", "HMCN1", "NPR3", "POSTN"],
        "FB": ["DCN", "ABCA8", "ABCA6", "ABCA10", "FBLN1", "COL15A1", "FBN1", "C7"],
        "L": ["SKAP1", "RIPOR2", "CD247", "IKZF1", "BCL11B", "SLFN12L", "ITGAL", "SAMD3", "CARD11", "CDC42SE2", "CCND3"],
        "MESO": ["C3", "SULF1", "AP000561.1", "PRG4", "GPM6A", "CDON", "DPP6", "CCDC80", "EZR", "FOS", "BNC1", "AC245041.2", "PRKD1", "CYSTM1", "TLL1", "WT1"],
        "MP": ["TBXAS1", "SLC9A9", "MRC1", "MS4A6A", "RBM47", "DOCK2", "MCTP1", "SYK", "MSR1", "ATP8B4", "F13A1", "CD74", "MS4A4E", "ADAP2"],
        "CM_cyto": ["TTN", "RYR2", "PAM", "TNNT2", "RABGAP1L", "PDLIM5", "MYL7", "MYH6"],
        "CM_nucl": ["RBM20", "TECRL", "MLIP", "CHRM2", "TRDN", "PALLD", "SGCD", "CMYA5", "MYOM2", "TBX5", "ESRRG", "LINC02248", "KCNJ3", "TACC2", "CORIN"],
    }
    
    print("    Creating cell type marker dotplots...")
    for cell_type, markers in marker_dict.items():
        available_markers = [gene for gene in markers if gene in adata_copy.var_names]
        
        if available_markers and len(available_markers) >= 3:  # At least 3 markers
            try:
                fig = sc.pl.dotplot(
                    adata_copy,
                    var_names=available_markers,
                    groupby='leiden',
                    dendrogram=True,
                    standard_scale='var',
                    figsize=(max(10, len(available_markers) * 0.4), 6),
                    cmap='Reds',
                    title=f'{cell_type} ({len(available_markers)} genes)',
                    show=False,
                    return_fig=True
                )
                if output_dir:
                    safe_name = cell_type.replace(" ", "_").replace("/", "_")
                    fig.savefig(os.path.join(output_dir, f'dotplot_{safe_name}.png'), 
                               dpi=150, bbox_inches='tight')
                plt.show()
                plt.close()
            except Exception as e:
                print(f"      Warning: Could not create dotplot for {cell_type}: {e}")
    
    # Hierarchical clustering
    print("    Performing hierarchical clustering...")
    for clst in np.unique(adata_copy.obs['leiden']):
        clst_size = sum(adata_copy.obs['leiden'] == clst)
        sc.tl.leiden(adata_copy, restrict_to=('leiden', [clst]), 
                    resolution=resolution_function(clst_size), key_added='leiden_R')
        adata_copy.obs['leiden'] = adata_copy.obs['leiden_R']
    
    # Calculate cluster-level scrublet scores
    print("    Calculating cluster-level doublet scores...")
    clst_meds = []
    for clst in np.unique(adata_copy.obs['leiden']):
        k = adata_copy.obs['leiden'] == clst
        clst_med = np.median(adata_copy.obs.loc[k, 'scrublet_score'])
        adata_copy.obs.loc[k, 'cluster_scrublet_score'] = clst_med
        clst_meds.append(clst_med)
    
    clst_meds = np.array(clst_meds)
    pvals, bh_pvals = test_outlier(clst_meds)
    
    for i, clst in enumerate(np.unique(adata_copy.obs['leiden'])):
        k = adata_copy.obs['leiden'] == clst
        adata_copy.obs.loc[k, 'pval'] = pvals[i]
        adata_copy.obs.loc[k, 'bh_pval'] = bh_pvals[i]
    
    # Create final UMAP plots
    print("    Creating doublet p-value UMAP plots...")
    fig2 = plt.figure(figsize=(16, 6), dpi=150)
    sc.pl.umap(adata_copy, color=['bh_pval','cluster_scrublet_score'],
               show=False, title=['Benjamini-Hochberg p-value', 'Cluster Scrublet Score'])
    plt.tight_layout()
    if output_dir:
        plt.savefig(os.path.join(output_dir, 'umap_pvals.png'), dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()
    
    # Transfer results back to original adata
    print("    Transferring results to original data...")
    adata.obs['scrublet_leiden'] = adata_copy.obs['leiden']
    adata.obs['scrublet_score'] = adata_copy.obs['scrublet_score']
    adata.obs['cluster_scrublet_score'] = adata_copy.obs['cluster_scrublet_score']
    adata.obs['doublet_pval'] = adata_copy.obs['pval']
    adata.obs['doublet_bh_pval'] = adata_copy.obs['bh_pval']
    
    # Identify doublets based on BH p-value
    adata.obs['predicted_doublet'] = adata.obs['doublet_bh_pval'] <= 0.05
    
    n_doublets = adata.obs['predicted_doublet'].sum()
    print(f"    Detected {n_doublets} doublets ({100*n_doublets/adata.n_obs:.1f}%) using BH p-value < 0.05")
    
    del adata_copy
    sc.settings.verbosity = old_verbosity
    
    return adata, scrub


def create_scrublet_plots(scrub, adata):
    """Create Scrublet diagnostic plots - removed doublet histogram"""
    plots = {}
    
    # Note: Doublet histogram removed as requested
    
    return plots


def save_html_report(plots_dict, output_path, sample_name, stats, output_dir):
    """Save all plots as a single HTML report with embedded images"""
    import base64
    from io import BytesIO
    from PIL import Image
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>QC Report - {sample_name}</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 20px;
                background-color: #f5f5f5;
            }}
            h1 {{
                color: #333;
                border-bottom: 2px solid #4CAF50;
                padding-bottom: 10px;
            }}
            h2 {{
                color: #555;
                margin-top: 30px;
                border-bottom: 1px solid #ddd;
                padding-bottom: 5px;
            }}
            .stats {{
                background-color: white;
                padding: 15px;
                border-radius: 5px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                margin: 20px 0;
            }}
            .stats table {{
                width: 100%;
                border-collapse: collapse;
            }}
            .stats td {{
                padding: 8px;
                border-bottom: 1px solid #eee;
            }}
            .stats td:first-child {{
                font-weight: bold;
                width: 40%;
            }}
            .plot {{
                background-color: white;
                padding: 15px;
                margin: 20px 0;
                border-radius: 5px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                text-align: center;
            }}
            .plot img {{
                max-width: 100%;
                height: auto;
            }}
            .warning {{
                background-color: #fff3cd;
                border-left: 4px solid #ffc107;
                padding: 10px;
                margin: 20px 0;
            }}
            .timestamp {{
                color: #888;
                font-size: 0.9em;
                margin-top: 20px;
            }}
        </style>
    </head>
    <body>
        <h1>Quality Control Report: {sample_name}</h1>
        <p class="timestamp">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        
        <div class="stats">
            <h2>Sample Statistics</h2>
            <table>
                <tr><td>Total cells (initial)</td><td>{stats['n_cells_initial']}</td></tr>
                <tr><td>After outlier filtering</td><td>{stats['n_cells_after_outliers']}</td></tr>
                <tr><td>After doublet removal</td><td>{stats['n_cells_final']}</td></tr>
                <tr><td>Total cells removed</td><td>{stats['n_cells_removed']} ({stats['pct_removed']:.1f}%)</td></tr>
                <tr><td>Outliers removed</td><td>{stats['n_outliers']} ({stats['pct_outliers']:.1f}%)</td></tr>
                <tr><td>MT outliers 3MAD (not removed)</td><td>{stats['n_mt_outliers']} ({stats['pct_mt_outliers']:.1f}%)</td></tr>
                <tr><td>Doublets removed</td><td>{stats['n_doublets']} ({stats['pct_doublets']:.1f}%)</td></tr>
                <tr><td>SoupX contamination</td><td>{stats['soupx_contamination']:.2%}</td></tr>
                <tr><td>Median genes per cell (final)</td><td>{stats['median_genes']:.0f}</td></tr>
                <tr><td>Median UMIs per cell (final)</td><td>{stats['median_umis']:.0f}</td></tr>
                <tr><td>Median MT% (final)</td><td>{stats['median_mt']:.2f}%</td></tr>
            </table>
        </div>
    """
    
    # Add warning if doublet rate is high
    if stats['pct_doublets'] > 15:
        html_content += f"""
        <div class="warning">
            <strong>⚠️ Warning:</strong> High doublet detection rate ({stats['pct_doublets']:.1f}%). 
            This may indicate sample quality issues or aggressive doublet detection parameters.
        </div>
        """
    
    # Function to convert plot to base64
    def fig_to_base64(fig):
        buf = BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
        buf.seek(0)
        img_str = base64.b64encode(buf.read()).decode()
        plt.close(fig)
        return img_str
    
    # Function to convert image file to base64
    def file_to_base64(filepath):
        try:
            with open(filepath, 'rb') as f:
                img_str = base64.b64encode(f.read()).decode()
            return img_str
        except Exception as e:
            print(f"Warning: Could not read {filepath}: {e}")
            return None
    
    # Add plots to HTML
    section_titles = {
        'before_qc': 'Before QC',
        'after_outliers': 'After Outlier Removal (Before SoupX & Doublets)',
        'after_soupx': 'After SoupX Correction',
        'scrublet': 'Doublet Detection',
        'final': 'Final (After All QC)'
    }
    
    for section, plots in plots_dict.items():
        if plots:
            html_content += f"<h2>{section_titles.get(section, section)}</h2>"
            for plot_name, fig in plots.items():
                img_str = fig_to_base64(fig)
                html_content += f"""
                <div class="plot">
                    <h3>{plot_name.replace('_', ' ').title()}</h3>
                    <img src="data:image/png;base64,{img_str}">
                </div>
                """
    
    # Add SoupX contamination plot
    soupx_plot_path = os.path.join(output_dir, 'soupx_contamination_plot.png')
    if os.path.exists(soupx_plot_path):
        html_content += f"<h2>SoupX Contamination Removal</h2>"
        img_str = file_to_base64(soupx_plot_path)
        if img_str:
            html_content += f"""
            <div class="plot">
                <h3>SoupX Contamination Plot</h3>
                <img src="data:image/png;base64,{img_str}">
            </div>
            """
    
    # Add UMAP leiden scrublet plot
    umap_leiden_path = os.path.join(output_dir, 'umap_leiden_scrublet.png')
    if os.path.exists(umap_leiden_path):
        html_content += f"<h2>Clustering and Doublet Scores</h2>"
        img_str = file_to_base64(umap_leiden_path)
        if img_str:
            html_content += f"""
            <div class="plot">
                <h3>UMAP - Leiden Clusters and Scrublet Scores</h3>
                <img src="data:image/png;base64,{img_str}">
            </div>
            """
    
    # Add cell type marker dotplots
    dotplot_files = sorted([f for f in os.listdir(output_dir) if f.startswith('dotplot_') and f.endswith('.png')])
    if dotplot_files:
        html_content += f"<h2>Cell Type Marker Expression</h2>"
        for dotplot_file in dotplot_files:
            dotplot_path = os.path.join(output_dir, dotplot_file)
            cell_type = dotplot_file.replace('dotplot_', '').replace('.png', '').replace('_', ' ')
            img_str = file_to_base64(dotplot_path)
            if img_str:
                html_content += f"""
                <div class="plot">
                    <h3>{cell_type}</h3>
                    <img src="data:image/png;base64,{img_str}">
                </div>
                """
    
    # Add UMAP p-values plot
    umap_pvals_path = os.path.join(output_dir, 'umap_pvals.png')
    if os.path.exists(umap_pvals_path):
        html_content += f"<h2>Doublet Detection P-values</h2>"
        img_str = file_to_base64(umap_pvals_path)
        if img_str:
            html_content += f"""
            <div class="plot">
                <h3>UMAP - Doublet P-values and Cluster Scores</h3>
                <img src="data:image/png;base64,{img_str}">
            </div>
            """
    
    html_content += """
    </body>
    </html>
    """
    
    with open(output_path, 'w') as f:
        f.write(html_content)
    
    print(f"  HTML report saved: {output_path}")


def process_sample(h5ad_path, raw_matrix_path, output_dir, sample_name):
    """Main processing function for one sample"""
    
    print(f"\n{'='*80}")
    print(f"Processing: {sample_name}")
    print(f"{'='*80}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    print("Loading data...")
    adata = sc.read_h5ad(h5ad_path)
    n_cells_initial = adata.n_obs
    
    # Calculate QC metrics
    print("Calculating QC metrics...")
    adata.var['MT'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['MT'], percent_top=[20], log1p=True, inplace=True
    )
    
    # Rename for consistency
    if 'pct_counts_MT' not in adata.obs.columns:
        adata.obs['pct_counts_MT'] = adata.obs['pct_counts_mt']
    
    # Create before QC plots
    print("Creating before-QC plots...")
    plots_before = create_qc_plots(adata, title_prefix="Before QC: ")
    
    # Detect outliers
    print("Detecting outliers...")
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    n_outliers = adata.obs.outlier.sum()
    print(f"  General outliers: {n_outliers} ({100*n_outliers/n_cells_initial:.1f}%)")
    
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_MT", 3) | (
        adata.obs["pct_counts_MT"] > 10
    )
    n_mt_outliers = adata.obs.mt_outlier.sum()
    print(f"  MT outliers: {n_mt_outliers} ({100*n_mt_outliers/n_cells_initial:.1f}%)")
    
    # Filter outliers
    print("\nFiltering outliers...")
    adata = adata[(~adata.obs.outlier) & (adata.obs.pct_counts_MT < 5)].copy()
    n_cells_after_outliers = adata.n_obs
    print(f"  Cells after outlier filtering: {n_cells_after_outliers}")
    
    # Create plots after outlier removal
    print("Creating after-outlier plots...")
    plots_after_outliers = create_qc_plots(adata, title_prefix="After Outliers: ")
    
    # Run SoupX - REQUIRED, do not continue if it fails
    adata = run_soupx(adata, raw_matrix_path, sample_name, output_dir)
    soupx_contamination = adata.uns.get('soupx_contamination', 0)
    
    # Create plots after SoupX
    print("Creating after-SoupX plots...")
    plots_after_soupx = create_qc_plots(adata, title_prefix="After SoupX: ")
    
    # Filter genes
    print(f"\nFiltering genes (min 20 cells)...")
    print(f"  Total number of genes before: {adata.n_vars}")
    sc.pp.filter_genes(adata, min_cells=20)
    print(f"  Number of genes after: {adata.n_vars}")
    
    # Run Scrublet with clustering
    adata, scrub = run_scrublet(adata, output_dir=output_dir)
    scrublet_plots = create_scrublet_plots(scrub, adata)
    
    # Remove doublets based on BH p-value
    print("\nRemoving doublets (BH p-value ≤ 0.05)...")
    n_doublets = adata.obs['predicted_doublet'].sum()
    adata_final = adata[~adata.obs['predicted_doublet']].copy()
    n_cells_final = adata_final.n_obs
    
    print(f"  Original cells (after outliers): {adata.n_obs}")
    print(f"  After removing doublets: {n_cells_final}")
    print(f"  Doublets removed: {n_doublets}")
    
    # Create final plots
    print("Creating final plots...")
    plots_final = create_qc_plots(adata_final, title_prefix="Final: ")
    
    # Collect statistics
    stats = {
        'n_cells_initial': n_cells_initial,
        'n_cells_after_outliers': n_cells_after_outliers,
        'n_cells_final': n_cells_final,
        'n_cells_removed': n_cells_initial - n_cells_final,
        'pct_removed': 100 * (n_cells_initial - n_cells_final) / n_cells_initial,
        'n_outliers': n_outliers,
        'pct_outliers': 100 * n_outliers / n_cells_initial,
        'n_mt_outliers': n_mt_outliers,
        'pct_mt_outliers': 100 * n_mt_outliers / n_cells_initial,
        'n_doublets': n_doublets,
        'pct_doublets': 100 * n_doublets / n_cells_after_outliers,
        'soupx_contamination': soupx_contamination,
        'median_genes': np.median(adata_final.obs['n_genes_by_counts']),
        'median_umis': np.median(adata_final.obs['total_counts']),
        'median_mt': np.median(adata_final.obs['pct_counts_MT'])
    }
    
    # Save outputs
    print("\nSaving outputs...")
    
    # Save filtered h5ad
    output_h5ad = os.path.join(output_dir, f"{sample_name}_QC_filtered.h5ad")
    adata_final.write_h5ad(output_h5ad)
    print(f"  Filtered data saved: {output_h5ad}")
    
    # Save HTML report
    all_plots = {
        'before_qc': plots_before,
        'after_outliers': plots_after_outliers,
        'after_soupx': plots_after_soupx,
        'scrublet': scrublet_plots,
        'final': plots_final
    }
    
    output_html = os.path.join(output_dir, f"{sample_name}_QC_report.html")
    save_html_report(all_plots, output_html, sample_name, stats, output_dir)
    
    print(f"\n✓ Sample {sample_name} processing complete!")
    print(f"{'='*80}\n")
    
    return adata_final, stats


# Main execution
if __name__ == "__main__":
    # Base directory
    base_dir = "/work/archive/public_studies/Hill/cellranger_runs/"
    output_base = "/work/archive/public_studies/Hill/qc_outputs/"
    
    # Get all sample directories
    all_items = os.listdir(base_dir)
    sample_dirs = [item for item in all_items 
                   if os.path.isdir(os.path.join(base_dir, item)) 
                   and not item.startswith('$')]
    
    print(f"Found {len(sample_dirs)} samples to process")
    
    # Process each sample
    all_stats = []
    failed_samples = []
    
    for i, sample_dir in enumerate(sample_dirs, 1):
        print(f"\n{'#'*80}")
        print(f"Sample {i}/{len(sample_dirs)}: {sample_dir}")
        print(f"{'#'*80}")
        
        try:
            # Construct paths
            h5ad_path = os.path.join(base_dir, sample_dir, sample_dir, 
                                     "outs", f"{sample_dir}_QClus.h5ad")
            raw_matrix_path = os.path.join(base_dir, sample_dir, sample_dir,
                                          "outs", "raw_feature_bc_matrix.h5")
            
            if not os.path.exists(h5ad_path):
                print(f"⚠️  Warning: File not found: {h5ad_path}")
                failed_samples.append({'sample': sample_dir, 'reason': 'h5ad not found'})
                continue
                
            if not os.path.exists(raw_matrix_path):
                print(f"⚠️  Warning: Raw matrix not found: {raw_matrix_path}")
                failed_samples.append({'sample': sample_dir, 'reason': 'raw matrix not found'})
                continue
            
            output_dir = os.path.join(output_base, sample_dir)
            
            # Process sample
            adata_filtered, stats = process_sample(h5ad_path, raw_matrix_path, 
                                                   output_dir, sample_dir)
            stats['sample_name'] = sample_dir
            all_stats.append(stats)
            
        except Exception as e:
            print(f"❌ ERROR processing {sample_dir}: {str(e)}")
            import traceback
            traceback.print_exc()
            failed_samples.append({'sample': sample_dir, 'reason': str(e)})
            continue
    
    # Save summary statistics
    if all_stats:
        summary_df = pd.DataFrame(all_stats)
        summary_path = os.path.join(output_base, "QC_summary.csv")
        summary_df.to_csv(summary_path, index=False)
        print(f"\n✓ Summary statistics saved: {summary_path}")
        
        # Print summary
        print(f"\n{'='*80}")
        print("PROCESSING SUMMARY")
        print(f"{'='*80}")
        print(f"Successfully processed: {len(all_stats)}/{len(sample_dirs)} samples")
        print(f"Failed: {len(failed_samples)}/{len(sample_dirs)} samples")
        
        if failed_samples:
            print("\nFailed samples:")
            for fail in failed_samples:
                print(f"  - {fail['sample']}: {fail['reason']}")
        
        print(f"\nOverall statistics:")
        print(f"  Mean cells removed: {summary_df['pct_removed'].mean():.1f}%")
        print(f"  Mean doublet rate: {summary_df['pct_doublets'].mean():.1f}%")
        print(f"  Mean SoupX contamination: {summary_df['soupx_contamination'].mean():.2%}")
    
    print(f"\n{'='*80}")
    print("All samples processed!")
    print(f"{'='*80}")
