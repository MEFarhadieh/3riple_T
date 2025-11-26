import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

sc.settings.set_figure_params(dpi=150, frameon=False, figsize=(8, 6))
sc.settings.verbosity = 1

adata = sc.read_h5ad('/work/archive/public_studies/merged_data/HB_integration_results/HB_integrated_data.h5ad')
adata

marker_dict = {
    "Nuclear genes (nucl_30)": ['MALAT1', 'NEAT1', 'FTX', 'FOXP1', 'RBMS3', 'ZBTB20', 'LRMDA', 'PBX1', 'ITPR2', 'AUTS2', 'TTC28', 'BNC2', 'EXOC4', 'RORA', 'PRKG1', 'ARID1B', 'PARD3B', 'GPHN', 'N4BP2L2', 'PKHD1L1', 'EXOC6B', 'FBXL7', 'MED13L', 'TBC1D5', 'IMMP2L', 'SYNE1', 'RERE', 'MBD5', 'EXT1', 'WWOX'],
    "Mitochondrial genes": ['MT-ND1', 'MT-ND2', 'MT-CO1', 'MT-CO2', 'MT-ATP8', 'MT-ATP6', 'MT-CO3', 'MT-ND3', 'MT-ND4L', 'MT-ND4', 'MT-ND5', 'MT-ND6', 'MT-CYB'],
    "VEC (Vascular Endothelial Cells)": ["VWF", "ERG", "ANO2", "PTPRB", "EGFL7", "PREX2", "ADGRL4", "FLT1", "CYYR1", "GRB10", "PPP1R16B", "DOCK9", "SHANK3", "PECAM1", "PLEKHG1", "EMCN"],
    "PER (Pericytes)": ["RGS5", "DACH1", "GUCY1A1", "ABCC9", "BGN", "NOTCH3", "PDGFRB", "FRMD3", "RNF152", "CCDC102B", "NGF"],
    "SMC (Smooth Muscle Cells)": ["MYH11", "KCNAB1", "NTRK3", "CHRM3", "ACTA2", "RGS6", "DGKG", "ITGA8", "TBX2", "LMOD1", "SDK1", "GPC6", "ANTXR1", "FLNA", "CLMN", "ATP10A", "MCAM", "TAGLN", "CCDC3"],
    "AD (Adipocytes)": ["PLIN4", "PLIN1", "PDE3B", "GPAM", "PTPRS", "PPARG", "MLXIPL", "MGST1", "AQP7", "SLC19A3", "FABP4", "TPRG1", "DIRC3", "LPL", "PNPLA2", "LIPE", "ADH1B", "ADIPOQ", "PRKAR2B", "CIDEA", "LINC00278", "PFKFB3", "LINC02237", "LIPE-AS1", "SVEP1"],
    "SC (Schwann Cells)": ["XKR4", "AC016766.1", "SLC35F1", "ZNF536", "NCAM2", "GPM6B", "KIRREL3", "SORCS1", "ST6GALNAC5", "PRKCA", "GINS3", "PMP22", "ALDH1A1", "IL1RAPL2", "DOCK5", "NKAIN3", "COL28A1", "RALGPS2", "PKN2-AS1", "KLHL29", "PTPRZ1"],
    "N (Neuronal)": ["CSMD1", "SYT1", "KCNIP4", "CNTNAP2", "DLGAP1", "PTPRD", "LRRTM4", "ATRNL1", "LRP1B", "CTNND2", "KCNQ5", "NRG3", "SNTG1", "GRIA2", "RIMS2", "CSMD3", "XIST", "KAZN", "DPP10", "HS6ST3", "OPCML"],
    "EEC (Endocardial Endothelial Cells)": ["PCDH7", "PCDH15", "LINC02147", "LINC02388", "MYRIP", "GMDS", "ADAMTSL1", "LEPR", "CALCRL", "CGNL1", "HMCN1", "NPR3", "POSTN"],
    "FB (Fibroblasts)": ["DCN", "ABCA8", "ABCA6", "ABCA10", "FBLN1", "COL15A1", "FBN1", "C7"],
    "L (Lymphocytes)": ["SKAP1", "RIPOR2", "CD247", "IKZF1", "BCL11B", "SLFN12L", "ITGAL", "SAMD3", "CARD11", "CDC42SE2", "CCND3"],
    "MESO (Mesothelial)": ["C3", "SULF1", "AP000561.1", "PRG4", "GPM6A", "CDON", "DPP6", "CCDC80", "EZR", "FOS", "BNC1", "AC245041.2", "PRKD1", "CYSTM1", "TLL1", "WT1"],
    "MP (Macrophages)": ["TBXAS1", "SLC9A9", "MRC1", "MS4A6A", "RBM47", "DOCK2", "MCTP1", "SYK", "MSR1", "ATP8B4", "F13A1", "CD74", "MS4A4E", "ADAP2"],
    "CM_cyto (Cardiomyocytes - Cytoplasmic)": ["TTN", "RYR2", "PAM", "TNNT2", "RABGAP1L", "PDLIM5", "MYL7", "MYH6"],
    "CM_nucl (Cardiomyocytes - Nuclear)": ["RBM20", "TECRL", "MLIP", "CHRM2", "TRDN", "PALLD", "SGCD", "CMYA5", "MYOM2", "TBX5", "ESRRG", "LINC02248", "KCNJ3", "TACC2", "CORIN", "DPY19L2", "WNK2", "MITF", "OBSCN", "FHOD3", "MYLK3", "DAPK2", "NEXN"],
    "T_B_Mast_MK (T/B/Mast/Megakaryocytes)": ['IL7R', 'CD247', 'SEL1L3', 'TPD52', 'CPAZ', 'KIT', 'NRGN', 'CD226']
}

for cell_type, markers in marker_dict.items():
    print(f"\n{'='*60}")
    print(f"Processing: {cell_type}")
    print(f"{'='*60}")
    
    available_markers = [gene for gene in markers if gene in adata.var_names]
    missing_markers = [gene for gene in markers if gene not in adata.var_names]

    print(f"Total markers: {len(markers)}")
    print(f"Available: {len(available_markers)}")
    print(f"Missing: {len(missing_markers)}")
    
    if missing_markers:
        print(f"Missing genes: {', '.join(missing_markers)}")
    
 
    if available_markers:
        print(f"\nDrawing dotplot for {len(available_markers)} genes...")
        
        #  dotplot
        sc.pl.dotplot(
            adata,
            var_names=available_markers,
            groupby='leiden_r0.8',
            dendrogram=True,  
            standard_scale='var',  
            figsize=(max(10, len(available_markers) * 0.4), 6),  
            cmap='Reds',  
            title=f'{cell_type}\n({len(available_markers)} genes)'
        )
        plt.tight_layout()
        plt.show()
    else:
        print(f"⚠️ WARNING: No genes available for {cell_type}!")
    
    print(f"{'='*60}\n")

print("\n✅ All dotplots completed!")

plt.tight_layout()
plt.show()

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity

# ========== 1st: Unknown Clusters ==========
print("=" * 60)
print("Step 1: Removing Unknown Clusters")
print("=" * 60)

unknown_clusters = ['14']
n_before = adata.n_obs

adata_filtered = adata[~adata.obs['leiden_r0.8'].isin(unknown_clusters)].copy()

adata_filtered.obs['leiden_r0.8'] = adata_filtered.obs['leiden_r0.8'].astype(str)

print(f"Removed {n_before - adata_filtered.n_obs:,} cells from unknown clusters")
print(f"Remaining: {adata_filtered.n_obs:,} cells\n")

# ========== 2nd: Sub-clustering Doublet Clusters ==========
print("=" * 60)
print("Step 2: Sub-clustering Doublet Clusters")
print("=" * 60)

doublet_clusters = ['18', '16', '17']

for cluster in doublet_clusters:
    print(f"\nProcessing doublet cluster {cluster}...")
    
    cluster_cells = adata_filtered[adata_filtered.obs['leiden_r0.8'] == cluster].copy()
    n_cells = cluster_cells.n_obs
    print(f"  Cells in cluster {cluster}: {n_cells}")
    
    if n_cells < 10:
        print(f"  Too few cells, skipping...")
        continue
    
    # Re-run neighbors and clustering for this cluster only
    sc.pp.neighbors(cluster_cells, n_neighbors=min(15, n_cells//2), n_pcs=30)
    sc.tl.leiden(cluster_cells, resolution=2.0, key_added='sub_cluster')
    
    n_subclusters = cluster_cells.obs['sub_cluster'].nunique()
    print(f"  Split into {n_subclusters} sub-clusters")
    
    # Markers for sub-cluster
    if n_subclusters > 1:
        # Differential expression between sub-clusters
        sc.tl.rank_genes_groups(cluster_cells, groupby='sub_cluster', method='wilcoxon')
        
        # Visulaize top genes for first 10 genes
        print(f"  Top marker genes per sub-cluster (showing first 10):")
        for subclust in sorted(cluster_cells.obs['sub_cluster'].unique())[:10]:
            markers = sc.get.rank_genes_groups_df(cluster_cells, group=subclust)
            top_genes = markers.head(5)['names'].tolist()
            print(f"    Sub {subclust}: {', '.join(top_genes)}")
    
    new_labels = cluster + '_' + cluster_cells.obs['sub_cluster'].astype(str)
    adata_filtered.obs.loc[cluster_cells.obs_names, 'leiden_r0.8'] = new_labels.values
    
    # Visualization
    sc.tl.umap(cluster_cells)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    sc.pl.umap(cluster_cells, color='sub_cluster', ax=axes[0], show=False,
               title=f'Cluster {cluster} - Sub-clusters', legend_loc='right margin',
               legend_fontsize=6)
    sc.pl.umap(cluster_cells, color='scrublet_score', ax=axes[1], show=False,
               title=f'Cluster {cluster} - Scrublet Score', cmap='RdYlBu_r')
    
    plt.tight_layout()
    plt.savefig(f'doublet_cluster_{cluster}_subclustering.png', dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()

print("\n✓ Doublet clusters sub-clustered")

# ========== 3rd: Analysis Bridge Clusters ==========
print("\n" + "=" * 60)
print("Step 3: Removing Bridge Cells from Bridge Clusters")
print("=" * 60)

bridge_clusters = ['19', '13', '12', '10', '9', '6']

#  Bridge cells identifier
def identify_bridge_cells(adata_cluster, density_threshold_percentile=10):
    """
    Identify bridge cells based on  UMAP density
    """
    
    # Calculate UMAP if not avail
    if 'X_umap' not in adata_cluster.obsm:
        sc.pp.neighbors(adata_cluster, n_neighbors=15, n_pcs=30)
        sc.tl.umap(adata_cluster)
    
    umap_coords = adata_cluster.obsm['X_umap']
    
    # Calculate density
    kde = KernelDensity(bandwidth=1.0, kernel='gaussian')
    kde.fit(umap_coords)
    log_density = kde.score_samples(umap_coords)
    density = np.exp(log_density)
    
    # Threshold (bottom percentile = low density = bridge)
    threshold = np.percentile(density, density_threshold_percentile)
    
    is_bridge = density < threshold
    
    return is_bridge, density

cells_removed_per_cluster = {}
bridge_cells_to_remove = []

for cluster in bridge_clusters:
    print(f"\nProcessing bridge cluster {cluster}...")
    
    cluster_cells = adata_filtered[adata_filtered.obs['leiden_r0.8'] == cluster].copy()
    n_cells = cluster_cells.n_obs
    print(f"  Cells in cluster {cluster}: {n_cells}")
    
    if n_cells < 10:
        print(f"  Too few cells, skipping...")
        continue
    
    # Identify bridge cells
    is_bridge, density = identify_bridge_cells(cluster_cells, density_threshold_percentile=10)
    
    # Combine with scrublet score
    high_scrublet = cluster_cells.obs['scrublet_score'] > 0.25
    
    # Bridge cells = low density OR high scrublet
    bridge_mask = is_bridge | high_scrublet
    
    n_bridge = bridge_mask.sum()
    pct_bridge = 100 * n_bridge / n_cells
    
    print(f"  Bridge cells identified: {n_bridge} ({pct_bridge:.1f}%)")
    print(f"    - Low density: {is_bridge.sum()}")
    print(f"    - High scrublet: {high_scrublet.sum()}")
    
    cells_removed_per_cluster[cluster] = n_bridge
    
    # Visualization - FIX: replace boolean to categorical with labels 
    cluster_cells.obs['is_bridge_cat'] = cluster_cells.obs.apply(
        lambda x: 'Bridge' if bridge_mask[x.name] else 'Core', axis=1
    ).astype('category')
    cluster_cells.obs['umap_density'] = density
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    sc.pl.umap(cluster_cells, color='umap_density', ax=axes[0], show=False,
               title=f'Cluster {cluster} - UMAP Density', cmap='viridis')
    
    sc.pl.umap(cluster_cells, color='scrublet_score', ax=axes[1], show=False,
               title=f'Cluster {cluster} - Scrublet Score', cmap='RdYlBu_r')
    
    # Categorical palette 
    sc.pl.umap(cluster_cells, color='is_bridge_cat', ax=axes[2], show=False,
               title=f'Cluster {cluster} - Bridge Cells', 
               palette={'Core': 'lightgray', 'Bridge': 'red'})
    
    plt.tight_layout()
    plt.savefig(f'bridge_cluster_{cluster}_analysis.png', dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()
    
    
    bridge_cell_names = cluster_cells.obs_names[bridge_mask].tolist()
    bridge_cells_to_remove.extend(bridge_cell_names)

print("\n✓ Bridge cells identified")
print("\nSummary of bridge cells per cluster:")
for cluster, n_removed in cells_removed_per_cluster.items():
    print(f"  Cluster {cluster}: {n_removed} cells")

# Remove bridge cells 
print(f"\nRemoving {len(bridge_cells_to_remove)} total bridge cells...")
adata_filtered = adata_filtered[~adata_filtered.obs_names.isin(bridge_cells_to_remove)].copy()

# ========== 4th: Re-clustering  ==========
print("\n" + "=" * 60)
print("Step 4: Final Re-clustering")
print("=" * 60)

print(f"Total cells remaining: {adata_filtered.n_obs:,}")
print(f"Total cells removed: {n_before - adata_filtered.n_obs:,} ({100*(n_before - adata_filtered.n_obs)/n_before:.2f}%)")

# Re-clustering
print("\nRe-computing neighbors, UMAP, and clustering...")
sc.pp.neighbors(adata_filtered, n_neighbors=30, n_pcs=20, use_rep='X_pca_harmony')
sc.tl.umap(adata_filtered, min_dist=0.3, spread=1.0)
sc.tl.leiden(adata_filtered, resolution=1.5, key_added='leiden_final')

print(f"Final number of clusters: {adata_filtered.obs['leiden_final'].nunique()}")

# ========== 5th: Visualization  ==========
print("\n" + "=" * 60)
print("Step 5: Final Visualization")
print("=" * 60)

adata_filtered.obs['leiden_r0.8'] = adata_filtered.obs['leiden_r0.8'].astype('category')

fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# Original
sc.pl.umap(adata, color='leiden_r0.8', ax=axes[0, 0], show=False,
           title=f'Original ({adata.obs["leiden_r0.8"].nunique()} clusters)',
           legend_loc='on data', legend_fontsize=6)

# After sub-clustering doublets
sc.pl.umap(adata_filtered, color='leiden_r0.8', ax=axes[0, 1], show=False,
           title=f'After sub-clustering & bridge removal',
           legend_loc='on data', legend_fontsize=5)

# Final
sc.pl.umap(adata_filtered, color='leiden_final', ax=axes[1, 0], show=False,
           title=f'Final ({adata_filtered.obs["leiden_final"].nunique()} clusters)',
           legend_loc='on data', legend_fontsize=6)

# Scrublet score
sc.pl.umap(adata_filtered, color='scrublet_score', ax=axes[1, 1], show=False,
           title='Scrublet Score (final)', cmap='RdYlBu_r', vmax=0.5)

plt.tight_layout()
plt.savefig('final_comparison.png', dpi=150, bbox_inches='tight')
plt.show()
plt.close()

# ========== 6th: Quality Check ==========
print("\n" + "=" * 60)
print("Step 6: Quality Check")
print("=" * 60)

# Scrublet distribution -  top 30 cluster
fig, axes = plt.subplots(1, 2, figsize=(18, 5))

# Before - top 30 clusters
import seaborn as sns
top_clusters_before = adata.obs['leiden_r0.8'].value_counts().head(30).index
data_before = adata.obs[adata.obs['leiden_r0.8'].isin(top_clusters_before)]
sns.violinplot(data=data_before, x='leiden_r0.8', y='scrublet_score', ax=axes[0])
axes[0].set_title('Scrublet Score - Original (top 30 clusters)')
axes[0].axhline(y=0.25, color='red', linestyle='--', label='Threshold')
axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=90, fontsize=6)
axes[0].legend()

# After - top 30 clusters
top_clusters_after = adata_filtered.obs['leiden_final'].value_counts().head(30).index
data_after = adata_filtered.obs[adata_filtered.obs['leiden_final'].isin(top_clusters_after)]
sns.violinplot(data=data_after, x='leiden_final', y='scrublet_score', ax=axes[1])
axes[1].set_title('Scrublet Score - Final (top 30 clusters)')
axes[1].axhline(y=0.25, color='red', linestyle='--', label='Threshold')
axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=90, fontsize=6)
axes[1].legend()

plt.tight_layout()
plt.savefig('scrublet_comparison.png', dpi=150, bbox_inches='tight')
plt.show()
plt.close()

# Statistics
print("\nFinal statistics:")
print(f"  Mean scrublet score: {adata_filtered.obs['scrublet_score'].mean():.3f}")
print(f"  Median scrublet score: {adata_filtered.obs['scrublet_score'].median():.3f}")
print(f"  % cells with score > 0.25: {(adata_filtered.obs['scrublet_score'] > 0.25).sum() / adata_filtered.n_obs * 100:.1f}%")

# Cluster size distribution
print("\nTop 20 largest clusters:")
cluster_sizes = adata_filtered.obs['leiden_final'].value_counts().head(20)
print(cluster_sizes)

# Cell typing
cluster_to_celltype = {
    '0': 'FB',
    '1': 'FB',
    '2': 'FB',
    '3': 'CM',
    '4': 'PER',
    '5': 'MP',
    '6': 'VEC',
    '7': 'MP',
    '8': 'CM',
    '9': 'AD', 
    '10': 'SMC',
    '11': 'FB',
    '12': 'TL',
    '13': 'CM',
    '14': 'EEC',
    '15': 'VEC',
    '16': 'VEC',
    '17': 'VEC',
    '18': 'SC',
    '19': 'MESO',
    '20': 'MP',
    '21': 'CM',
    '22': 'FB',
    '23': 'MP',
    '24': 'BL',
    '25': 'CM',
    '26': 'FB',
    '27': 'MAST',
    '28': 'CM',
   
}


adata_filtered.obs['cell_type'] = adata_filtered.obs['res_final'].map(cluster_to_celltype)

unmapped_clusters = adata_filtered.obs[adata_filtered.obs['cell_type'].isna()]['res_final'].unique()
if len(unmapped_clusters) > 0:
    print(f"Note: These clusters are unmapped: {unmapped_clusters}")

print("\nNumber of cells in each cell type:")
print(adata_filtered.obs['cell_type'].value_counts().sort_index())

sc.pl.umap(adata_filtered, color='cell_type', legend_loc='on data', 
           title='Cell Type Annotation', frameon=False, 
           legend_fontsize=10, legend_fontoutline=2, size=3)

sc.pl.umap(adata_filtered, color='res_final', legend_loc='right margin', frameon=False, size=3)

# ========== 7th: Save  ==========
print("\n" + "=" * 60)
print("Step 7: Saving Results")
print("=" * 60)

adata_filtered.write('/work/archive/public_studies/merged_data/HB_integration_results/HB_adata_cleaned_final.h5ad')
print("✓ Saved: /work/archive/public_studies/merged_data/HB_integration_results/HB_adata_cleaned_final.h5ad")
