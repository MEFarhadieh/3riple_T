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

# Define datasets with their paths and study names
DATASETS = {
    'Hill': '/work/archive/public_studies/Hill/qc_outputs/',
    'Bondareva_SR': '/work/archive/public_studies/Bondareva_SR/qc_outputs/',
    'Bondareva_AF': '/work/archive/public_studies/Bondareva_AF/qc_outputs/',
    'Bondareva_pAF': '/work/archive/public_studies/Bondareva_pAF/qc_outputs/'
}

# Map to study names (Bondareva variants all mapped to 'Bondareva')
STUDY_MAPPING = {
    'Hill': 'Hill',
    'Bondareva_SR': 'Bondareva',
    'Bondareva_AF': 'Bondareva',
    'Bondareva_pAF': 'Bondareva',
}

def load_sample(h5ad_path, sample_name, dataset_name, study_name):
    """Load a single sample and add metadata"""
    print(f"  Loading {sample_name}...")
    
    try:
        adata = sc.read_h5ad(h5ad_path)
        
        # Add metadata
        adata.obs['sample_id'] = sample_name
        adata.obs['dataset'] = dataset_name
        adata.obs['study'] = study_name
        
        # Store sample info in uns
        if 'sample_info' not in adata.uns:
            adata.uns['sample_info'] = {}
        
        adata.uns['sample_info']['sample_id'] = sample_name
        adata.uns['sample_info']['dataset'] = dataset_name
        adata.uns['sample_info']['study'] = study_name
        
        print(f"    Loaded: {adata.n_obs} cells, {adata.n_vars} genes")
        
        return adata
        
    except Exception as e:
        print(f"    ERROR loading {sample_name}: {e}")
        return None


def collect_all_samples():
    """Collect all samples from all datasets"""
    print("="*80)
    print("COLLECTING ALL SAMPLES")
    print("="*80)
    
    all_samples = []
    dataset_summary = []
    
    for dataset_name, base_path in DATASETS.items():
        print(f"\nDataset: {dataset_name}")
        print(f"Path: {base_path}")
        
        study_name = STUDY_MAPPING[dataset_name]
        
        if not os.path.exists(base_path):
            print(f"  WARNING: Path does not exist!")
            continue
        
        # Get all sample directories
        sample_dirs = [d for d in os.listdir(base_path) 
                      if os.path.isdir(os.path.join(base_path, d))]
        
        print(f"  Found {len(sample_dirs)} sample directories")
        
        n_loaded = 0
        for sample_dir in sample_dirs:
            h5ad_path = os.path.join(base_path, sample_dir, 
                                    f"{sample_dir}_QC_filtered.h5ad")
            
            if not os.path.exists(h5ad_path):
                print(f"    WARNING: {sample_dir} - h5ad not found")
                continue
            
            adata = load_sample(h5ad_path, sample_dir, dataset_name, study_name)
            
            if adata is not None:
                all_samples.append(adata)
                n_loaded += 1
                
                # Track for summary
                dataset_summary.append({
                    'sample_id': sample_dir,
                    'dataset': dataset_name,
                    'study': study_name,
                    'n_cells': adata.n_obs,
                    'n_genes': adata.n_vars
                })
        
        print(f"  Successfully loaded: {n_loaded}/{len(sample_dirs)} samples")
    
    return all_samples, pd.DataFrame(dataset_summary)


def merge_samples(adata_list):
    """Merge all samples into one AnnData object"""
    print("\n" + "="*80)
    print("MERGING SAMPLES")
    print("="*80)
    
    print(f"\nTotal samples to merge: {len(adata_list)}")
    
    # Check total cells
    total_cells = sum([adata.n_obs for adata in adata_list])
    print(f"Total cells: {total_cells:,}")
    
    # Merge with outer join (keep all genes)
    print("\nMerging datasets (outer join)...")
    adata_merged = sc.concat(
        adata_list,
        join='outer',
        label='batch',
        keys=[adata.obs['sample_id'].iloc[0] for adata in adata_list],
        index_unique='_'
    )
    
    print(f"\nMerged dataset:")
    print(f"  Cells: {adata_merged.n_obs:,}")
    print(f"  Genes: {adata_merged.n_vars:,}")
    
    # Fill NaN values with 0 (for genes not present in all samples)
    if hasattr(adata_merged.X, 'data'):
        # Sparse matrix - already handled
        pass
    else:
        # Dense matrix
        adata_merged.X = np.nan_to_num(adata_merged.X, 0)
    
    print("\nSamples per study:")
    print(adata_merged.obs['study'].value_counts())
    
    print("\nSamples per dataset:")
    print(adata_merged.obs['dataset'].value_counts())
    
    return adata_merged


def normalize_and_preprocess(adata):
    """Normalize, find HVGs, and compute PCA - PRESERVE ALL LAYERS"""
    print("\n" + "="*80)
    print("NORMALIZATION AND PREPROCESSING")
    print("="*80)
    
    # Store raw counts if not already stored
    if 'counts' not in adata.layers:
        print("\nStoring raw counts in layer 'counts'...")
        adata.layers['counts'] = adata.X.copy()
    
    # If soupX_counts exists, preserve it
    if 'soupX_counts' in adata.layers:
        print("Preserving soupX_counts layer...")
    
    print(f"\nInitial data:")
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    print(f"  Layers: {list(adata.layers.keys())}")
    
    # Normalize to 10,000 counts per cell
    print("\nNormalizing to 10,000 counts per cell...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Log transform
    print("Log1p transformation...")
    sc.pp.log1p(adata)
    
    # Store normalized in layer
    print("Storing log-normalized data in layer 'log1p_norm'...")
    adata.layers['log1p_norm'] = adata.X.copy()
    
    # Highly variable genes - batch-aware
    print("\nFinding highly variable genes (batch-aware)...")
    sc.pp.highly_variable_genes(
        adata,
        batch_key='sample_id',
        n_top_genes=4000,
        subset=False,  # Don't subset, just flag
        flavor='seurat_v3'
    )
    
    n_hvg = adata.var['highly_variable'].sum()
    print(f"  Identified {n_hvg} highly variable genes")
    
    # Scale data (on HVGs only for PCA)
    print("\nScaling data...")
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    
    # PCA
    print("\nComputing PCA (50 components)...")
    sc.tl.pca(adata_hvg, n_comps=50, svd_solver='arpack')
    
    # Transfer PCA back to original object
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.varm['PCs'] = np.zeros((adata.n_vars, 50))
    adata.varm['PCs'][adata.var['highly_variable'], :] = adata_hvg.varm['PCs']
    adata.uns['pca'] = adata_hvg.uns['pca']
    
    del adata_hvg
    
    print("\nPreprocessing complete!")
    print(f"  Final layers: {list(adata.layers.keys())}")
    
    return adata


def create_qc_plots(adata, output_dir):
    """Create QC plots for merged dataset"""
    print("\n" + "="*80)
    print("CREATING QC PLOTS")
    print("="*80)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Plot 1: Cells per study
    fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
    study_counts = adata.obs['study'].value_counts()
    ax.bar(study_counts.index, study_counts.values)
    ax.set_xlabel('Study')
    ax.set_ylabel('Number of Cells')
    ax.set_title('Cells per Study')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cells_per_study.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: cells_per_study.png")
    
    # Plot 2: Cells per dataset
    fig, ax = plt.subplots(figsize=(12, 6), dpi=150)
    dataset_counts = adata.obs['dataset'].value_counts()
    ax.bar(range(len(dataset_counts)), dataset_counts.values)
    ax.set_xticks(range(len(dataset_counts)))
    ax.set_xticklabels(dataset_counts.index, rotation=45, ha='right')
    ax.set_xlabel('Dataset')
    ax.set_ylabel('Number of Cells')
    ax.set_title('Cells per Dataset')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cells_per_dataset.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: cells_per_dataset.png")
    
    # Plot 3: Cells per sample
    fig, ax = plt.subplots(figsize=(20, 6), dpi=150)
    sample_counts = adata.obs['sample_id'].value_counts()
    ax.bar(range(len(sample_counts)), sample_counts.values)
    ax.set_xlabel('Sample')
    ax.set_ylabel('Number of Cells')
    ax.set_title(f'Cells per Sample (n={len(sample_counts)})')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cells_per_sample.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: cells_per_sample.png")
    
    # Plot 4: QC metrics by study
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=150)
    
    # Total counts
    sns.violinplot(data=adata.obs, x='study', y='total_counts', ax=axes[0,0])
    axes[0,0].set_title('Total Counts per Study')
    axes[0,0].set_xlabel('')
    
    # Number of genes
    sns.violinplot(data=adata.obs, x='study', y='n_genes_by_counts', ax=axes[0,1])
    axes[0,1].set_title('Number of Genes per Study')
    axes[0,1].set_xlabel('')
    
    # MT percentage
    if 'pct_counts_MT' in adata.obs.columns:
        sns.violinplot(data=adata.obs, x='study', y='pct_counts_MT', ax=axes[1,0])
        axes[1,0].set_title('MT% per Study')
        axes[1,0].set_xlabel('Study')
    
    # Cells per study (bar)
    study_counts = adata.obs['study'].value_counts()
    axes[1,1].bar(range(len(study_counts)), study_counts.values)
    axes[1,1].set_xticks(range(len(study_counts)))
    axes[1,1].set_xticklabels(study_counts.index, rotation=45, ha='right')
    axes[1,1].set_ylabel('Number of Cells')
    axes[1,1].set_title('Cells per Study')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'qc_metrics_by_study.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: qc_metrics_by_study.png")
    
    # Plot 5: PCA variance
    fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
    variance_ratio = adata.uns['pca']['variance_ratio']
    ax.plot(range(1, len(variance_ratio) + 1), 
            np.cumsum(variance_ratio) * 100, 'o-')
    ax.set_xlabel('PC')
    ax.set_ylabel('Cumulative Variance Explained (%)')
    ax.set_title('PCA Cumulative Variance')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pca_variance.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: pca_variance.png")
    
    # Plot 6: PCA by study
    fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
    studies = adata.obs['study'].unique()
    colors = plt.cm.tab10(range(len(studies)))
    
    for i, study in enumerate(studies):
        mask = adata.obs['study'] == study
        ax.scatter(adata.obsm['X_pca'][mask, 0],
                  adata.obsm['X_pca'][mask, 1],
                  c=[colors[i]], label=study, s=1, alpha=0.5)
    
    ax.set_xlabel(f'PC1 ({variance_ratio[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2 ({variance_ratio[1]*100:.1f}%)')
    ax.set_title('PCA colored by Study')
    ax.legend(markerscale=5)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pca_by_study.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: pca_by_study.png")


def save_summary_stats(adata, summary_df, output_dir):
    """Save summary statistics"""
    print("\n" + "="*80)
    print("SAVING SUMMARY STATISTICS")
    print("="*80)
    
    # Sample summary
    summary_df.to_csv(os.path.join(output_dir, 'samples_summary.csv'), index=False)
    print(f"  Saved: samples_summary.csv")
    
    # Study summary
    study_summary = adata.obs.groupby('study').agg({
        'sample_id': 'nunique',
        'total_counts': ['median', 'mean'],
        'n_genes_by_counts': ['median', 'mean']
    }).round(2)
    study_summary.columns = ['n_samples', 'median_counts', 'mean_counts', 
                             'median_genes', 'mean_genes']
    study_summary['n_cells'] = adata.obs['study'].value_counts()
    study_summary.to_csv(os.path.join(output_dir, 'study_summary.csv'))
    print(f"  Saved: study_summary.csv")
    
    # Dataset summary
    dataset_summary = adata.obs.groupby('dataset').agg({
        'sample_id': 'nunique',
        'total_counts': ['median', 'mean'],
        'n_genes_by_counts': ['median', 'mean']
    }).round(2)
    dataset_summary.columns = ['n_samples', 'median_counts', 'mean_counts',
                               'median_genes', 'mean_genes']
    dataset_summary['n_cells'] = adata.obs['dataset'].value_counts()
    dataset_summary.to_csv(os.path.join(output_dir, 'dataset_summary.csv'))
    print(f"  Saved: dataset_summary.csv")
    
    # Overall stats
    overall_stats = {
        'total_cells': adata.n_obs,
        'total_genes': adata.n_vars,
        'n_studies': adata.obs['study'].nunique(),
        'n_datasets': adata.obs['dataset'].nunique(),
        'n_samples': adata.obs['sample_id'].nunique(),
        'n_hvgs': adata.var['highly_variable'].sum(),
        'median_counts': adata.obs['total_counts'].median(),
        'median_genes': adata.obs['n_genes_by_counts'].median(),
        'layers': list(adata.layers.keys())
    }
    
    with open(os.path.join(output_dir, 'overall_stats.txt'), 'w') as f:
        f.write("OVERALL STATISTICS\n")
        f.write("="*50 + "\n\n")
        for key, value in overall_stats.items():
            f.write(f"{key}: {value}\n")
    print(f"  Saved: overall_stats.txt")
    
    # Print summary
    print("\n" + "="*80)
    print("OVERALL SUMMARY")
    print("="*80)
    for key, value in overall_stats.items():
        print(f"  {key}: {value}")


def main():
    """Main execution"""
    print("\n" + "#"*80)
    print("# MERGE AND PREPROCESS ALL DATASETS")
    print("#"*80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # Output directory
    output_base = "/work/archive/public_studies/merged_data/"
    os.makedirs(output_base, exist_ok=True)
    
    # Step 1: Collect all samples
    adata_list, summary_df = collect_all_samples()
    
    if len(adata_list) == 0:
        print("\nERROR: No samples loaded!")
        return
    
    print(f"\nTotal samples loaded: {len(adata_list)}")
    
    # Step 2: Merge
    adata_merged = merge_samples(adata_list)
    
    # Clear memory
    del adata_list
    
    # Step 3: Normalize and preprocess
    adata_merged = normalize_and_preprocess(adata_merged)
    
    # Step 4: Create QC plots
    create_qc_plots(adata_merged, output_base)
    
    # Step 5: Save summary stats
    save_summary_stats(adata_merged, summary_df, output_base)
    
    # Step 6: Save merged and preprocessed data
    print("\n" + "="*80)
    print("SAVING MERGED DATA")
    print("="*80)
    
    output_h5ad = os.path.join(output_base, 'HB_merged_preprocessed.h5ad')
    print(f"\nSaving to: {output_h5ad}")
    adata_merged.write_h5ad(output_h5ad)
    
    # Calculate file size
    file_size_gb = os.path.getsize(output_h5ad) / (1024**3)
    print(f"File size: {file_size_gb:.2f} GB")
    
    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"\nOutput directory: {output_base}")
    print(f"Merged data: {output_h5ad}")
    print(f"  - Cells: {adata_merged.n_obs:,}")
    print(f"  - Genes: {adata_merged.n_vars:,}")
    print(f"  - Layers: {list(adata_merged.layers.keys())}")
    

if __name__ == "__main__":
    main()
