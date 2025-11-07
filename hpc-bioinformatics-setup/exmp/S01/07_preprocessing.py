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


def test_outlier_improved(x, upper_mad_only=True):
    """
    Improved outlier detection with better handling of edge cases
    """
    
    # Check minimum data points
    if len(x) < 3:
        print(f"    Warning: Only {len(x)} clusters - insufficient for statistical testing")
        return np.ones(len(x)), np.ones(len(x))
    
    med = np.median(x)
    
    # Calculate MAD
    if upper_mad_only:
        upper_values = x[x > med]
        if len(upper_values) < 2:
            # Not enough upper values, use all
            mad = np.median(np.abs(x - med)) * 1.4826
        else:
            mad = np.median(upper_values - med) * 1.4826
    else:
        mad = np.median(np.abs(x - med)) * 1.4826
    
    # Handle edge case: MAD too small (all values very similar)
    if mad < 0.01:  # threshold for "too similar"
        print(f"    Warning: MAD={mad:.4f} is very small. Cluster scores are too similar for statistical testing.")
        print(f"    Score range: [{x.min():.3f}, {x.max():.3f}]")
        # Return high p-values (not significant)
        return np.ones(len(x)), np.ones(len(x))
    
    # Calculate p-values using normal distribution
    pvals = 1 - scipy.stats.norm.cdf(x, loc=med, scale=mad)
    
    # Check if all p-values are too similar
    if np.std(pvals) < 0.01:
        print(f"    Warning: All p-values are very similar (std={np.std(pvals):.4f})")
    
    # Apply BH correction
    try:
        bh_pvals = multipletests(pvals, method='fdr_bh')[1]
    except:
        bh_pvals = pvals.copy()
    
    return pvals, bh_pvals


def detect_doublet_clusters_hybrid(adata_copy):
    """
    Hybrid approach: combines statistical testing with absolute thresholds
    This is more robust for different sample qualities
    """
    
    print("    Calculating cluster-level doublet scores...")
    
    clst_meds = []
    cluster_ids = []
    cluster_sizes = []
    
    for clst in np.unique(adata_copy.obs['leiden']):
        k = adata_copy.obs['leiden'] == clst
        clst_med = np.median(adata_copy.obs.loc[k, 'scrublet_score'])
        clst_size = k.sum()
        
        adata_copy.obs.loc[k, 'cluster_scrublet_score'] = clst_med
        clst_meds.append(clst_med)
        cluster_ids.append(clst)
        cluster_sizes.append(clst_size)
    
    clst_meds = np.array(clst_meds)
    cluster_sizes = np.array(cluster_sizes)
    n_clusters = len(clst_meds)
    
    print(f"    Number of clusters: {n_clusters}")
    print(f"    Cluster score range: [{clst_meds.min():.3f}, {clst_meds.max():.3f}]")
    print(f"    Cluster score median: {np.median(clst_meds):.3f}")
    print(f"    Cluster score 75th percentile: {np.percentile(clst_meds, 75):.3f}")
    
    # Method 1: Statistical testing (MAD-based)
    pvals, bh_pvals = test_outlier_improved(clst_meds, upper_mad_only=True)
    
    # Method 2: Percentile-based (relative threshold)
    # Clusters in top 25% of scores
    score_75th = np.percentile(clst_meds, 75)
    percentile_mask = clst_meds > score_75th
    
    # Method 3: Absolute threshold
    # Clusters with very high scrublet scores
    absolute_threshold = 0.3  # You can adjust this
    absolute_mask = clst_meds > absolute_threshold
    
    # Method 4: Z-score based (robust to outliers)
    z_scores = (clst_meds - np.median(clst_meds)) / (np.std(clst_meds) + 1e-10)
    zscore_mask = z_scores > 2  # 2 standard deviations above median
    
    # Combine methods
    # Strategy: A cluster is doublet if it meets ANY of these criteria:
    # 1. Statistical test passes (BH p-value < 0.05)
    # 2. In top 25% AND above absolute threshold
    # 3. Very high Z-score (> 2)
    
    doublet_statistical = bh_pvals < 0.05
    doublet_combined = percentile_mask & absolute_mask
    doublet_extreme = zscore_mask
    
    final_doublet_mask = doublet_statistical | doublet_combined | doublet_extreme
    
    # Store all results in adata
    for i, clst in enumerate(cluster_ids):
        k = adata_copy.obs['leiden'] == clst
        adata_copy.obs.loc[k, 'pval'] = pvals[i]
        adata_copy.obs.loc[k, 'bh_pval'] = bh_pvals[i]
        adata_copy.obs.loc[k, 'z_score'] = z_scores[i]
        adata_copy.obs.loc[k, 'is_doublet_statistical'] = doublet_statistical[i]
        adata_copy.obs.loc[k, 'is_doublet_combined'] = doublet_combined[i]
        adata_copy.obs.loc[k, 'is_doublet_extreme'] = doublet_extreme[i]
        adata_copy.obs.loc[k, 'is_doublet_final'] = final_doublet_mask[i]
    
    # Report
    n_doublet_stat = doublet_statistical.sum()
    n_doublet_comb = doublet_combined.sum()
    n_doublet_extr = doublet_extreme.sum()
    n_doublet_final = final_doublet_mask.sum()
    
    print(f"\n    Doublet cluster detection results:")
    print(f"      Method 1 - Statistical (BH p<0.05):     {n_doublet_stat}/{n_clusters} clusters")
    print(f"      Method 2 - Combined threshold:          {n_doublet_comb}/{n_clusters} clusters")
    print(f"      Method 3 - Extreme Z-score (>2):        {n_doublet_extr}/{n_clusters} clusters")
    print(f"      FINAL - Any method:                     {n_doublet_final}/{n_clusters} clusters")
    
    if n_doublet_final > 0:
        doublet_clusters = [cluster_ids[i] for i in range(n_clusters) if final_doublet_mask[i]]
        doublet_scores = [clst_meds[i] for i in range(n_clusters) if final_doublet_mask[i]]
        print(f"\n    Doublet clusters: {doublet_clusters}")
        print(f"    Their scores: {[f'{s:.3f}' for s in doublet_scores]}")
    
    return adata_copy, pvals, bh_pvals


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
    
    # Use hybrid doublet detection method
    adata_copy, pvals, bh_pvals = detect_doublet_clusters_hybrid(adata_copy)
    
    # Create final UMAP plots - UPDATED TO SAVE SEPARATELY
    print("    Creating doublet p-value UMAP plots...")
    
    # Combined plot (existing)
    fig2 = plt.figure(figsize=(16, 6), dpi=150)
    sc.pl.umap(adata_copy, color=['bh_pval','cluster_scrublet_score'],
               show=False, title=['Benjamini-Hochberg p-value', 'Cluster Scrublet Score'])
    plt.tight_layout()
    if output_dir:
        plt.savefig(os.path.join(output_dir, 'umap_pvals.png'), dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()
    
    # Individual UMAP plots for HTML report
    if output_dir: 
        # Z-score
        fig_z = plt.figure(figsize=(8, 6), dpi=150)
        sc.pl.umap(adata_copy, color='z_score', show=False,
                   title='Z-score', cmap='RdYlBu_r')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'umap_z_score.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
        # Final doublet prediction - CONVERT BOOLEAN TO LABELED CATEGORICAL
        # FIX: Convert boolean to categorical with proper labels
        adata_copy.obs['doublet_status'] = adata_copy.obs['is_doublet_final'].map({
            True: 'Doublet',
            False: 'Singlet'
        }).astype('category')
        
        fig_final = plt.figure(figsize=(8, 6), dpi=150)
        sc.pl.umap(adata_copy, color='doublet_status', show=False,
                   title='Final Doublet Prediction', 
                   palette={'Singlet': 'lightgray', 'Doublet': 'red'})
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'umap_doublet_final.png'), dpi=150, bbox_inches='tight')
        plt.close()
    
    # Transfer results back to original adata - SAFE VERSION
    print("    Transferring results to original data...")
    adata.obs['scrublet_leiden'] = adata_copy.obs['leiden'].values
    adata.obs['scrublet_score'] = adata_copy.obs['scrublet_score'].values
    adata.obs['cluster_scrublet_score'] = adata_copy.obs['cluster_scrublet_score'].values
    adata.obs['doublet_pval'] = adata_copy.obs['pval'].values
    adata.obs['doublet_bh_pval'] = adata_copy.obs['bh_pval'].values
    adata.obs['doublet_z_score'] = adata_copy.obs['z_score'].values
    adata.obs['is_doublet_statistical'] = adata_copy.obs['is_doublet_statistical'].values
    adata.obs['is_doublet_combined'] = adata_copy.obs['is_doublet_combined'].values
    adata.obs['is_doublet_extreme'] = adata_copy.obs['is_doublet_extreme'].values
    
    # CRITICAL FIX: Convert boolean to actual boolean values
    adata.obs['predicted_doublet'] = adata_copy.obs['is_doublet_final'].values.astype(bool)
    
    n_doublets = int(adata.obs['predicted_doublet'].sum())
    print(f"    Final: Detected {n_doublets} doublet cells ({100*n_doublets/adata.n_obs:.1f}%)")
    
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
    
    print(f"  DEBUG: Starting HTML report generation...")
    print(f"  DEBUG: Stats type: {type(stats)}")
    print(f"  DEBUG: Stats keys: {stats.keys() if isinstance(stats, dict) else 'Not a dict'}")
    
    # Ensure all stats are proper types and extract safely
    try:
        n_cells_initial = int(stats.get('n_cells_initial', 0))
        n_cells_after_outliers = int(stats.get('n_cells_after_outliers', 0))
        n_cells_final = int(stats.get('n_cells_final', 0))
        n_cells_removed = int(stats.get('n_cells_removed', 0))
        pct_removed = float(stats.get('pct_removed', 0))
        n_outliers = int(stats.get('n_outliers', 0))
        pct_outliers = float(stats.get('pct_outliers', 0))
        n_mt_outliers = int(stats.get('n_mt_outliers', 0))
        pct_mt_outliers = float(stats.get('pct_mt_outliers', 0))
        n_doublets = int(stats.get('n_doublets', 0))
        pct_doublets = float(stats.get('pct_doublets', 0))
        n_doublets_statistical = int(stats.get('n_doublets_statistical', 0))
        pct_doublets_statistical = float(stats.get('pct_doublets_statistical', 0))
        n_doublets_combined = int(stats.get('n_doublets_combined', 0))
        pct_doublets_combined = float(stats.get('pct_doublets_combined', 0))
        n_doublets_extreme = int(stats.get('n_doublets_extreme', 0))
        pct_doublets_extreme = float(stats.get('pct_doublets_extreme', 0))
        soupx_contamination = float(stats.get('soupx_contamination', 0))
        median_genes = float(stats.get('median_genes', 0))
        median_umis = float(stats.get('median_umis', 0))
        median_mt = float(stats.get('median_mt', 0))
        
        print(f"  DEBUG: All stats extracted successfully")
        
    except Exception as e:
        print(f"  ERROR extracting stats: {e}")
        print(f"  DEBUG: Problematic stats: {stats}")
        raise
    
    try:
        html_content = f"""<!DOCTYPE html>
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
        h3 {{
            color: #666;
            margin-top: 20px;
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
        .doublet-methods {{
            background-color: #e3f2fd;
            padding: 15px;
            border-radius: 5px;
            border-left: 4px solid #2196F3;
            margin: 20px 0;
        }}
        .doublet-methods table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 10px;
        }}
        .doublet-methods th {{
            background-color: #2196F3;
            color: white;
            padding: 10px;
            text-align: left;
        }}
        .doublet-methods td {{
            padding: 8px;
            border-bottom: 1px solid #ddd;
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
        .plot-grid {{
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 20px;
            margin: 20px 0;
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
            <tr><td>Total cells (initial)</td><td>{n_cells_initial}</td></tr>
            <tr><td>After outlier filtering</td><td>{n_cells_after_outliers}</td></tr>
            <tr><td>After doublet removal</td><td>{n_cells_final}</td></tr>
            <tr><td>Total cells removed</td><td>{n_cells_removed} ({pct_removed:.1f}%)</td></tr>
            <tr><td>Outliers removed</td><td>{n_outliers} ({pct_outliers:.1f}%)</td></tr>
            <tr><td>MT outliers 3MAD (not removed)</td><td>{n_mt_outliers} ({pct_mt_outliers:.1f}%)</td></tr>
            <tr><td>Doublets removed (total)</td><td>{n_doublets} ({pct_doublets:.1f}%)</td></tr>
            <tr><td>SoupX contamination</td><td>{soupx_contamination:.2%}</td></tr>
            <tr><td>Median genes per cell (final)</td><td>{median_genes:.0f}</td></tr>
            <tr><td>Median UMIs per cell (final)</td><td>{median_umis:.0f}</td></tr>
            <tr><td>Median MT% (final)</td><td>{median_mt:.2f}%</td></tr>
        </table>
    </div>
    
    <div class="doublet-methods">
        <h2>📊 Doublet Detection Methods (Hybrid Approach)</h2>
        <p>Our doublet detection uses three complementary methods. A cell is flagged as doublet if it passes <strong>any</strong> of these criteria:</p>
        <table>
            <tr>
                <th>Method</th>
                <th>Description</th>
                <th>Cells Detected</th>
                <th>Percentage</th>
            </tr>
            <tr>
                <td><strong>Method 1: Statistical</strong></td>
                <td>MAD-based outlier detection with Benjamini-Hochberg correction (p &lt; 0.05)</td>
                <td>{n_doublets_statistical}</td>
                <td>{pct_doublets_statistical:.2f}%</td>
            </tr>
            <tr>
                <td><strong>Method 2: Combined Threshold</strong></td>
                <td>Clusters in top 25% of scrublet scores AND above 0.3 threshold</td>
                <td>{n_doublets_combined}</td>
                <td>{pct_doublets_combined:.2f}%</td>
            </tr>
            <tr>
                <td><strong>Method 3: Extreme Z-score</strong></td>
                <td>Clusters with Z-score &gt; 2 (2 standard deviations above median)</td>
                <td>{n_doublets_extreme}</td>
                <td>{pct_doublets_extreme:.2f}%</td>
            </tr>
            <tr style="background-color: #f0f0f0; font-weight: bold;">
                <td colspan="2"><strong>TOTAL (Union of all methods)</strong></td>
                <td>{n_doublets}</td>
                <td>{pct_doublets:.2f}%</td>
            </tr>
        </table>
        <p style="margin-top: 10px; font-size: 0.9em; color: #555;">
            <strong>Note:</strong> The total is the final correct number and prcentage, while Cells Detected column is boolean, 0 is No and 1 is Yes.
        </p>
    </div>
"""
        
        print(f"  DEBUG: HTML header created successfully")
        
    except Exception as e:
        print(f"  ERROR creating HTML header: {e}")
        import traceback
        traceback.print_exc()
        raise
    
    # Add warning if doublet rate is high
    if pct_doublets > 15:
        html_content += f"""
    <div class="warning">
        <strong>⚠️ Warning:</strong> High doublet detection rate ({pct_doublets:.1f}%). 
        This may indicate sample quality issues or aggressive doublet detection parameters.
    </div>
"""
    
    # Function to convert plot to base64
    def fig_to_base64(fig):
        try:
            buf = BytesIO()
            fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
            buf.seek(0)
            img_str = base64.b64encode(buf.read()).decode('utf-8')
            plt.close(fig)
            return img_str
        except Exception as e:
            print(f"    Warning: Could not convert figure to base64: {e}")
            return None
    
    # Function to convert image file to base64
    def file_to_base64(filepath):
        try:
            with open(filepath, 'rb') as f:
                img_str = base64.b64encode(f.read()).decode('utf-8')
            return img_str
        except Exception as e:
            print(f"    Warning: Could not read {filepath}: {e}")
            return None
    
    # Add plots to HTML
    print(f"  DEBUG: Adding matplotlib plots...")
    section_titles = {
        'before_qc': 'Before QC',
        'after_outliers': 'After Outlier Removal (Before SoupX & Doublets)',
        'after_soupx': 'After SoupX Correction',
        'scrublet': 'Doublet Detection',
        'final': 'Final (After All QC)'
    }
    
    for section, plots in plots_dict.items():
        if plots:
            html_content += f"<h2>{section_titles.get(section, section)}</h2>\n"
            for plot_name, fig in plots.items():
                img_str = fig_to_base64(fig)
                if img_str:
                    html_content += f"""
    <div class="plot">
        <h3>{plot_name.replace('_', ' ').title()}</h3>
        <img src="data:image/png;base64,{img_str}">
    </div>
"""
    
    # Add SoupX contamination plot
    print(f"  DEBUG: Adding SoupX plot...")
    soupx_plot_path = os.path.join(output_dir, 'soupx_contamination_plot.png')
    if os.path.exists(soupx_plot_path):
        html_content += "<h2>SoupX Contamination Removal</h2>\n"
        img_str = file_to_base64(soupx_plot_path)
        if img_str:
            html_content += f"""
    <div class="plot">
        <h3>SoupX Contamination Plot</h3>
        <img src="data:image/png;base64,{img_str}">
    </div>
"""
    
    # Add UMAP leiden scrublet plot
    print(f"  DEBUG: Adding UMAP leiden plot...")
    umap_leiden_path = os.path.join(output_dir, 'umap_leiden_scrublet.png')
    if os.path.exists(umap_leiden_path):
        html_content += "<h2>Clustering and Doublet Scores</h2>\n"
        img_str = file_to_base64(umap_leiden_path)
        if img_str:
            html_content += f"""
    <div class="plot">
        <h3>UMAP - Leiden Clusters and Scrublet Scores</h3>
        <img src="data:image/png;base64,{img_str}">
    </div>
"""
    
    # Add cell type marker dotplots
    print(f"  DEBUG: Adding dotplots...")
    try:
        dotplot_files = sorted([f for f in os.listdir(output_dir) if f.startswith('dotplot_') and f.endswith('.png')])
        if dotplot_files:
            html_content += "<h2>Cell Type Marker Expression</h2>\n"
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
    except Exception as e:
        print(f"    Warning: Could not add dotplots: {e}")
    
    # Add detailed doublet detection UMAP plots
    print(f"  DEBUG: Adding detailed UMAP plots...")
    html_content += "<h2>Doublet Detection - Detailed UMAP Visualizations</h2>\n"
    
    # Grid of 4 UMAP plots
    umap_files = [
        ('umap_bh_pval.png', 'Benjamini-Hochberg P-value'),
        ('umap_cluster_score.png', 'Cluster Scrublet Score'),
        ('umap_z_score.png', 'Z-score'),
        ('umap_doublet_final.png', 'Final Doublet Prediction')
    ]
    
    html_content += '<div class="plot-grid">\n'
    for umap_file, title in umap_files:
        umap_path = os.path.join(output_dir, umap_file)
        if os.path.exists(umap_path):
            img_str = file_to_base64(umap_path)
            if img_str:
                html_content += f"""
    <div class="plot">
        <h3>{title}</h3>
        <img src="data:image/png;base64,{img_str}">
    </div>
"""
    html_content += '</div>\n'
    
    # Also add the combined plot
    umap_pvals_path = os.path.join(output_dir, 'umap_pvals.png')
    if os.path.exists(umap_pvals_path):
        img_str = file_to_base64(umap_pvals_path)
        if img_str:
            html_content += f"""
    <div class="plot">
        <h3>Combined View: P-values and Cluster Scores</h3>
        <img src="data:image/png;base64,{img_str}">
    </div>
"""
    
    html_content += """
</body>
</html>
"""
    
    print(f"  DEBUG: Writing HTML file...")
    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"  ✓ HTML report saved: {output_path}")
    except Exception as e:
        print(f"  ERROR writing HTML file: {e}")
        import traceback
        traceback.print_exc()
        raise

def process_sample(h5ad_path, raw_matrix_path, output_dir, sample_name):
    """Main processing function for one sample"""
    
    print(f"\n{'='*80}")
    print(f"Processing: {sample_name}")
    print(f"{'='*80}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    try:
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
        n_outliers = int(adata.obs.outlier.sum())
        print(f"  General outliers: {n_outliers} ({100*n_outliers/n_cells_initial:.1f}%)")
        
        adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_MT", 3) | (
            adata.obs["pct_counts_MT"] > 10
        )
        n_mt_outliers = int(adata.obs.mt_outlier.sum())
        print(f"  MT outliers: {n_mt_outliers} ({100*n_mt_outliers/n_cells_initial:.1f}%)")
        
        # Filter outliers - SAFE VERSION
        print("\nFiltering outliers...")
        # Create boolean mask safely
        keep_mask = (~adata.obs['outlier'].values) & (adata.obs['pct_counts_MT'].values < 5)
        adata = adata[keep_mask].copy()
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
        
        # Get detailed doublet statistics from different methods
        n_doublets_statistical = int(adata.obs.get('is_doublet_statistical', False).sum()) if 'is_doublet_statistical' in adata.obs.columns else 0
        n_doublets_combined = int(adata.obs.get('is_doublet_combined', False).sum()) if 'is_doublet_combined' in adata.obs.columns else 0
        n_doublets_extreme = int(adata.obs.get('is_doublet_extreme', False).sum()) if 'is_doublet_extreme' in adata.obs.columns else 0
        
        # Remove doublets based on hybrid method - SAFE VERSION
        print("\nRemoving doublets (hybrid method)...")
        # Convert to plain boolean array
        doublet_mask = adata.obs['predicted_doublet'].values.astype(bool)
        keep_mask = ~doublet_mask
        
        n_doublets = int(doublet_mask.sum())
        adata_final = adata[keep_mask].copy()
        n_cells_final = adata_final.n_obs
        
        print(f"  Original cells (after outliers): {adata.n_obs}")
        print(f"  After removing doublets: {n_cells_final}")
        print(f"  Doublets removed: {n_doublets}")
        
        # Create final plots
        print("Creating final plots...")
        plots_final = create_qc_plots(adata_final, title_prefix="Final: ")
        
        # Collect statistics - UPDATED WITH DETAILED DOUBLET INFO
        stats = {
            'n_cells_initial': int(n_cells_initial),
            'n_cells_after_outliers': int(n_cells_after_outliers),
            'n_cells_final': int(n_cells_final),
            'n_cells_removed': int(n_cells_initial - n_cells_final),
            'pct_removed': float(100 * (n_cells_initial - n_cells_final) / n_cells_initial),
            'n_outliers': int(n_outliers),
            'pct_outliers': float(100 * n_outliers / n_cells_initial),
            'n_mt_outliers': int(n_mt_outliers),
            'pct_mt_outliers': float(100 * n_mt_outliers / n_cells_initial),
            'n_doublets': int(n_doublets),
            'pct_doublets': float(100 * n_doublets / n_cells_after_outliers) if n_cells_after_outliers > 0 else 0,
            # Detailed doublet detection stats
            'n_doublets_statistical': int(n_doublets_statistical),
            'pct_doublets_statistical': float(100 * n_doublets_statistical / n_cells_after_outliers) if n_cells_after_outliers > 0 else 0,
            'n_doublets_combined': int(n_doublets_combined),
            'pct_doublets_combined': float(100 * n_doublets_combined / n_cells_after_outliers) if n_cells_after_outliers > 0 else 0,
            'n_doublets_extreme': int(n_doublets_extreme),
            'pct_doublets_extreme': float(100 * n_doublets_extreme / n_cells_after_outliers) if n_cells_after_outliers > 0 else 0,
            'soupx_contamination': float(soupx_contamination),
            'median_genes': float(np.median(adata_final.obs['n_genes_by_counts'])),
            'median_umis': float(np.median(adata_final.obs['total_counts'])),
            'median_mt': float(np.median(adata_final.obs['pct_counts_MT']))
        }
        
        # Save outputs
        print("\nSaving outputs...")
        
        # FIRST: Save HTML report (before h5ad which might fail)
        print("  Creating HTML report...")
        all_plots = {
            'before_qc': plots_before,
            'after_outliers': plots_after_outliers,
            'after_soupx': plots_after_soupx,
            'scrublet': scrublet_plots,
            'final': plots_final
        }
        
        output_html = os.path.join(output_dir, f"{sample_name}_QC_report.html")
        try:
            save_html_report(all_plots, output_html, sample_name, stats, output_dir)
        except Exception as e:
            print(f"  ERROR creating HTML report: {e}")
            import traceback
            traceback.print_exc()
        
        # THEN: Save filtered h5ad
        print("  Saving h5ad file...")
        
        # Fix boolean columns for h5ad compatibility
        technical_cols = ['is_doublet_statistical', 'is_doublet_combined', 
                         'is_doublet_extreme', 'outlier', 'mt_outlier',
                         'doublet_status', 'is_doublet_final_cat', 'is_doublet_final_int']
        
        for col in technical_cols:
            if col in adata_final.obs.columns:
                adata_final.obs = adata_final.obs.drop(columns=[col])
        
        # Convert remaining boolean columns to categorical
        for col in adata_final.obs.columns:
            if adata_final.obs[col].dtype == bool or adata_final.obs[col].dtype == 'bool':
                adata_final.obs[col] = adata_final.obs[col].map({True: 'True', False: 'False'}).astype('category')
        
        output_h5ad = os.path.join(output_dir, f"{sample_name}_QC_filtered.h5ad")
        try:
            adata_final.write_h5ad(output_h5ad)
            print(f"  ✓ Filtered data saved: {output_h5ad}")
        except Exception as e:
            print(f"  ERROR saving h5ad: {e}")
            import traceback
            traceback.print_exc()
        
    except Exception as e:
        print(f"\n❌ ERROR in sample processing: {str(e)}")
        import traceback
        traceback.print_exc()
        raise


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
