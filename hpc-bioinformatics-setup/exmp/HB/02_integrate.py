import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Set parameters
sc.settings.set_figure_params(dpi=150, frameon=False, figsize=(10, 8))
sc.settings.verbosity = 1

# Paths
merged_path = "/work/archive/public_studies/merged_data/HB_merged_preprocessed.h5ad"
output_dir = "/work/archive/public_studies/merged_data/HB_integration_results/"

os.makedirs(output_dir, exist_ok=True)

print("="*80)
print("INTEGRATION AND UMAP VISUALIZATION PIPELINE")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# ============================================================================
# 1. LOAD MERGED DATA
# ============================================================================
print("\n" + "="*80)
print("1. LOADING MERGED DATA")
print("="*80)

adata = sc.read_h5ad(merged_path)

print(f"\nLoaded data:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")
print(f"  Layers: {list(adata.layers.keys())}")

# Check batch information
print(f"\nBatch information:")
print(f"  Studies: {adata.obs['study'].nunique()}")
print(f"  Datasets: {adata.obs['dataset'].nunique()}")
print(f"  Samples: {adata.obs['sample_id'].nunique()}")

print("\nStudy distribution:")
print(adata.obs['study'].value_counts())

print("\nSample distribution:")
print(adata.obs['sample_id'].value_counts())

# ============================================================================
# 2. QUALITY CHECK
# ============================================================================
print("\n" + "="*80)
print("2. QUALITY METRICS SUMMARY")
print("="*80)

print(f"\nGlobal metrics:")
print(f"  Mean genes per cell: {adata.obs['n_genes_by_counts'].mean():.0f}")
print(f"  Median genes per cell: {adata.obs['n_genes_by_counts'].median():.0f}")
print(f"  Mean UMIs per cell: {adata.obs['total_counts'].mean():.0f}")
print(f"  Median UMIs per cell: {adata.obs['total_counts'].median():.0f}")
print(f"  HVGs: {adata.var['highly_variable'].sum()}")

# ============================================================================
# 3. PREPARE DATA FOR INTEGRATION
# ============================================================================
print("\n" + "="*80)
print("3. PREPARING DATA FOR INTEGRATION")
print("="*80)

# Use log1p_norm layer if available
if 'log1p_norm' in adata.layers:
    print("\nUsing log1p_norm layer...")
    adata.X = adata.layers['log1p_norm'].copy()
else:
    print("\nLog1p_norm layer not found, using current X...")

# Subset to HVGs for integration
print(f"\nSubsetting to {adata.var['highly_variable'].sum()} HVGs...")
adata_hvg = adata[:, adata.var['highly_variable']].copy()

print(f"HVG subset:")
print(f"  Cells: {adata_hvg.n_obs:,}")
print(f"  Genes: {adata_hvg.n_vars:,}")

# Scale data
print("\nScaling data...")
sc.pp.scale(adata_hvg, max_value=10)

# ============================================================================
# 4. PCA
# ============================================================================
print("\n" + "="*80)
print("4. COMPUTING PCA")
print("="*80)

print("\nRunning PCA (50 components)...")
sc.tl.pca(adata_hvg, n_comps=50, svd_solver='arpack')

# Plot PCA variance
print("\nCreating PCA variance plot...")
fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
variance_ratio = adata_hvg.uns['pca']['variance_ratio']
cumsum_var = np.cumsum(variance_ratio) * 100

ax.plot(range(1, 51), cumsum_var, 'o-', linewidth=2, markersize=4)
ax.axhline(y=50, color='r', linestyle='--', alpha=0.5, label='50% variance')
ax.axhline(y=80, color='orange', linestyle='--', alpha=0.5, label='80% variance')
ax.set_xlabel('Principal Component')
ax.set_ylabel('Cumulative Variance Explained (%)')
ax.set_title('PCA Cumulative Variance')
ax.grid(True, alpha=0.3)
ax.legend()
plt.tight_layout()
plt.savefig(f"{output_dir}/01_pca_variance_ratio.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: 01_pca_variance_ratio.png")

# Determine number of PCs to use
# Check if variance reaches 80%
if np.any(cumsum_var >= 80):
    n_pcs = np.where(cumsum_var >= 80)[0][0] + 1
    print(f"\nUsing {n_pcs} PCs (explains {cumsum_var[n_pcs-1]:.1f}% variance)")
else:
    # If doesn't reach 80%, use max variance point or 30 PCs
    n_pcs = 30
    print(f"\nVariance doesn't reach 80% with 50 PCs")
    print(f"Maximum variance explained: {cumsum_var[-1]:.1f}%")
    print(f"Using {n_pcs} PCs (explains {cumsum_var[n_pcs-1]:.1f}% variance)")

n_pcs = min(n_pcs, 30)  # Cap at 30 for efficiency

# ============================================================================
# 5. UMAP BEFORE INTEGRATION
# ============================================================================
print("\n" + "="*80)
print("5. UMAP BEFORE INTEGRATION")
print("="*80)

print("\nComputing neighbors...")
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=n_pcs)

print("Computing UMAP...")
sc.tl.umap(adata_hvg)

# Save UMAP before integration
print("\nCreating UMAP plots (before integration)...")

# By sample_id and study
fig, axes = plt.subplots(1, 2, figsize=(16, 6), dpi=150)
sc.pl.umap(adata_hvg, color='sample_id', ax=axes[0], show=False,
           title='Before Integration - by Sample ID', legend_loc='right margin',
           legend_fontsize=6, size=3)
sc.pl.umap(adata_hvg, color='study', ax=axes[1], show=False,
           title='Before Integration - by Study', legend_loc='on data',
           legend_fontsize=8, size=3)
plt.tight_layout()
plt.savefig(f"{output_dir}/02_umap_before_integration.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: 02_umap_before_integration.png")

# By QC metrics
fig, axes = plt.subplots(1, 3, figsize=(18, 5), dpi=150)
sc.pl.umap(adata_hvg, color='n_genes_by_counts', ax=axes[0], show=False,
           title='n_genes', cmap='viridis', size=3)
sc.pl.umap(adata_hvg, color='total_counts', ax=axes[1], show=False,
           title='total_counts', cmap='viridis', size=3)
if 'pct_counts_MT' in adata_hvg.obs.columns:
    sc.pl.umap(adata_hvg, color='pct_counts_MT', ax=axes[2], show=False,
               title='MT%', cmap='viridis', size=3)
plt.tight_layout()
plt.savefig(f"{output_dir}/03_umap_before_qc_metrics.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: 03_umap_before_qc_metrics.png")

# ============================================================================
# 6. HARMONY INTEGRATION
# ============================================================================
print("\n" + "="*80)
print("6. RUNNING HARMONY INTEGRATION")
print("="*80)

print("\nIntegrating by 'sample_id'...")
import scanpy.external as sce

sce.pp.harmony_integrate(
    adata_hvg,
    key='sample_id',
    basis='X_pca',
    adjusted_basis='X_pca_harmony',
    max_iter_harmony=20
)

print("  Harmony integration complete!")

# ============================================================================
# 7. UMAP AFTER INTEGRATION
# ============================================================================
print("\n" + "="*80)
print("7. UMAP AFTER INTEGRATION")
print("="*80)

print("\nComputing neighbors on integrated data...")
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=n_pcs, use_rep='X_pca_harmony')

print("Computing UMAP...")
sc.tl.umap(adata_hvg)

# Save UMAP after integration
print("\nCreating UMAP plots (after integration)...")

# By sample_id and study
fig, axes = plt.subplots(1, 2, figsize=(16, 6), dpi=150)
sc.pl.umap(adata_hvg, color='sample_id', ax=axes[0], show=False,
           title='After Integration - by Sample ID', legend_loc='right margin',
           legend_fontsize=6, size=3)
sc.pl.umap(adata_hvg, color='study', ax=axes[1], show=False,
           title='After Integration - by Study', legend_loc='on data',
           legend_fontsize=10, size=3)
plt.tight_layout()
plt.savefig(f"{output_dir}/04_umap_after_integration.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: 04_umap_after_integration.png")

# By QC metrics
fig, axes = plt.subplots(1, 3, figsize=(18, 5), dpi=150)
sc.pl.umap(adata_hvg, color='n_genes_by_counts', ax=axes[0], show=False,
           title='n_genes', cmap='viridis', size=3)
sc.pl.umap(adata_hvg, color='total_counts', ax=axes[1], show=False,
           title='total_counts', cmap='viridis', size=3)
if 'pct_counts_MT' in adata_hvg.obs.columns:
    sc.pl.umap(adata_hvg, color='pct_counts_MT', ax=axes[2], show=False,
               title='MT%', cmap='viridis', size=3)
plt.tight_layout()
plt.savefig(f"{output_dir}/05_umap_after_qc_metrics.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: 05_umap_after_qc_metrics.png")

# ============================================================================
# 8. CLUSTERING
# ============================================================================
print("\n" + "="*80)
print("8. CLUSTERING ON INTEGRATED DATA")
print("="*80)

resolutions = [0.8, 1, 1.5, 2]
print(f"\nTesting resolutions: {resolutions}")

for res in resolutions:
    sc.tl.leiden(adata_hvg, resolution=res, key_added=f'leiden_r{res}')
    n_clusters = adata_hvg.obs[f'leiden_r{res}'].nunique()
    print(f"  Resolution {res}: {n_clusters} clusters")

# Plot clustering results
print("\nCreating clustering plots...")
fig, axes = plt.subplots(2, 3, figsize=(18, 12), dpi=150)
axes = axes.flatten()

for i, res in enumerate(resolutions):
    sc.pl.umap(adata_hvg, color=f'leiden_r{res}', ax=axes[i], show=False,
               title=f'Leiden (resolution={res})',
               legend_loc='on data', legend_fontsize=6, size=3)

plt.tight_layout()
plt.savefig(f"{output_dir}/06_umap_clustering_resolutions.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: 06_umap_clustering_resolutions.png")

# ============================================================================
# 9. CLUSTER COMPOSITION ANALYSIS
# ============================================================================
print("\n" + "="*80)
print("9. CLUSTER COMPOSITION ANALYSIS")
print("="*80)

leiden_key = 'leiden_r0.8'
print(f"\nUsing {leiden_key} for composition analysis...")

# Sample composition per cluster
cluster_sample_comp = pd.crosstab(
    adata_hvg.obs[leiden_key],
    adata_hvg.obs['sample_id'],
    normalize='index'
) * 100

fig, ax = plt.subplots(figsize=(14, 8), dpi=150)
cluster_sample_comp.plot(kind='bar', stacked=True, ax=ax,
                         colormap='tab20', width=0.8)
ax.set_xlabel('Cluster', fontsize=12)
ax.set_ylabel('Percentage (%)', fontsize=12)
ax.set_title('Sample Composition per Cluster', fontsize=14)
ax.legend(title='Sample ID', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=6)
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.tight_layout()
plt.savefig(f"{output_dir}/07_cluster_sample_composition.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: 07_cluster_sample_composition.png")

# Study composition per cluster
cluster_study_comp = pd.crosstab(
    adata_hvg.obs[leiden_key],
    adata_hvg.obs['study'],
    normalize='index'
) * 100

fig, ax = plt.subplots(figsize=(14, 8), dpi=150)
cluster_study_comp.plot(kind='bar', stacked=True, ax=ax,
                        colormap='Set3', width=0.8)
ax.set_xlabel('Cluster', fontsize=12)
ax.set_ylabel('Percentage (%)', fontsize=12)
ax.set_title('Study Composition per Cluster', fontsize=14)
ax.legend(title='Study', bbox_to_anchor=(1.05, 1), loc='upper left')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.tight_layout()
plt.savefig(f"{output_dir}/08_cluster_study_composition.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: 08_cluster_study_composition.png")

# Dataset composition per cluster
cluster_dataset_comp = pd.crosstab(
    adata_hvg.obs[leiden_key],
    adata_hvg.obs['dataset'],
    normalize='index'
) * 100

fig, ax = plt.subplots(figsize=(14, 8), dpi=150)
cluster_dataset_comp.plot(kind='bar', stacked=True, ax=ax,
                          colormap='tab20c', width=0.8)
ax.set_xlabel('Cluster', fontsize=12)
ax.set_ylabel('Percentage (%)', fontsize=12)
ax.set_title('Dataset Composition per Cluster', fontsize=14)
ax.legend(title='Dataset', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.tight_layout()
plt.savefig(f"{output_dir}/09_cluster_dataset_composition.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: 09_cluster_dataset_composition.png")

# Cluster sizes
cluster_sizes = adata_hvg.obs[leiden_key].value_counts().sort_index()
fig, ax = plt.subplots(figsize=(12, 6), dpi=150)
ax.bar(cluster_sizes.index.astype(str), cluster_sizes.values)
ax.set_xlabel('Cluster')
ax.set_ylabel('Number of Cells')
ax.set_title(f'Cluster Sizes ({leiden_key})')
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig(f"{output_dir}/10_cluster_sizes.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: 10_cluster_sizes.png")

# ============================================================================
# 10. SAVE INTEGRATED DATA
# ============================================================================
print("\n" + "="*80)
print("10. SAVING INTEGRATED DATA")
print("="*80)

# Transfer integrated embeddings back to full adata
print("\nTransferring integrated embeddings to full dataset...")
adata.obsm['X_pca_harmony'] = np.zeros((adata.n_obs, n_pcs))
adata.obsm['X_pca_harmony'] = adata_hvg.obsm['X_pca_harmony'][:, :n_pcs]

adata.obsm['X_umap'] = adata_hvg.obsm['X_umap'].copy()

# Transfer clustering
for res in resolutions:
    adata.obs[f'leiden_r{res}'] = adata_hvg.obs[f'leiden_r{res}'].copy()

# Save
output_path = f"{output_dir}/HB_integrated_data.h5ad"
print(f"\nSaving to: {output_path}")
adata.write_h5ad(output_path)

file_size_gb = os.path.getsize(output_path) / (1024**3)
print(f"File size: {file_size_gb:.2f} GB")

# ============================================================================
# 11. SAVE SUMMARY STATISTICS
# ============================================================================
print("\n" + "="*80)
print("11. SAVING SUMMARY STATISTICS")
print("="*80)

# Integration summary
summary_stats = {
    'Total cells': f"{adata.n_obs:,}",
    'Total genes': f"{adata.n_vars:,}",
    'HVG genes used': f"{adata_hvg.n_vars:,}",
    'Studies': adata.obs['study'].nunique(),
    'Datasets': adata.obs['dataset'].nunique(),
    'Samples': adata.obs['sample_id'].nunique(),
    'PCs used': n_pcs,
    'Clusters (r=0.8)': adata_hvg.obs['leiden_r0.8'].nunique(),
    'Integration method': 'Harmony',
    'Batch key': 'sample_id',
    'Mean genes/cell': f"{adata.obs['n_genes_by_counts'].mean():.0f}",
    'Median genes/cell': f"{adata.obs['n_genes_by_counts'].median():.0f}",
    'Mean UMIs/cell': f"{adata.obs['total_counts'].mean():.0f}",
    'Date': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
}

summary_df = pd.DataFrame(list(summary_stats.items()),
                         columns=['Metric', 'Value'])
summary_df.to_csv(f"{output_dir}/integration_summary.csv", index=False)
print(f"  Saved: integration_summary.csv")

# Cluster composition tables
cluster_sample_comp.to_csv(f"{output_dir}/cluster_sample_composition.csv")
cluster_study_comp.to_csv(f"{output_dir}/cluster_study_composition.csv")
cluster_dataset_comp.to_csv(f"{output_dir}/cluster_dataset_composition.csv")
cluster_sizes.to_csv(f"{output_dir}/cluster_sizes.csv", header=['n_cells'])
print(f"  Saved: composition tables")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print("\n" + "="*80)
print("INTEGRATION COMPLETE!")
print("="*80)
print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

print("\n📊 Summary:")
for key, value in summary_stats.items():
    print(f"  {key}: {value}")

print(f"\n📁 Output directory: {output_dir}")
print("\n📄 Generated files:")
print("  1. HB_integrated_data.h5ad")
print("  2. 01_pca_variance_ratio.png")
print("  3. 02_umap_before_integration.png")
print("  4. 03_umap_before_qc_metrics.png")
print("  5. 04_umap_after_integration.png")
print("  6. 05_umap_after_qc_metrics.png")
print("  7. 06_umap_clustering_resolutions.png")
print("  8. 07_cluster_sample_composition.png")
print("  9. 08_cluster_study_composition.png")
print("  10. 09_cluster_dataset_composition.png")
print("  11. 10_cluster_sizes.png")
print("  12. integration_summary.csv")
print("  13. cluster_*_composition.csv")
print("  14. cluster_sizes.csv")

print("\n" + "="*80)
print("✅ ALL DONE!")
print("="*80)
