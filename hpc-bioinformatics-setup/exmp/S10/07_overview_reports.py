#!/usr/bin/env python3
"""
snRNA-seq QC Report Generator
Generates HTML QC reports for samples processed with QClus
Usage: python qc_reports.py <base_directory>
Example: python qc_reports.py /mnt/archive/farhadie/public_studies/Bondareva_AF/cellranger_runs
"""

import os
import sys
import glob
import base64
from io import BytesIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from datetime import datetime

# Settings
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)
plt.rcParams['figure.figsize'] = (10, 6)

def fig_to_base64(fig):
    """Convert matplotlib figure to base64 for embedding in HTML"""
    # Handle both Figure and Axes objects
    if hasattr(fig, 'figure'):
        fig = fig.figure
    elif isinstance(fig, list):
        # Some scanpy functions return list of axes
        fig = fig[0].figure if fig else plt.gcf()
    
    buf = BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=100)
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return f"data:image/png;base64,{img_base64}"

def determine_sex(adata):
    """Determine sex based on XIST expression"""
    if 'XIST' not in adata.var_names:
        return "Unknown", None
    
    passed_cells = adata.obs['qclus'] == 'passed'
    if passed_cells.sum() == 0:
        return "Unknown", None
    
    xist_expr = adata[passed_cells, 'XIST'].X.toarray().flatten()
    mean_xist = np.mean(xist_expr)
    
    sex = "Female" if mean_xist > 0.1 else "Male"
    return sex, xist_expr

def plot_xist_violin(adata):
    """XIST expression violin plot"""
    if 'XIST' not in adata.var_names:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, 'XIST gene not found', ha='center', va='center', fontsize=14)
        ax.axis('off')
        return fig
    
    desired_order = ['passed', 'scrublet filter', 'outlier filter', 'clustering filter', 'initial filter']
    qclus_order = [cat for cat in desired_order if cat in adata.obs['qclus'].unique()]
    
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.violin(adata, 'XIST', groupby='qclus', order=qclus_order, ax=ax, show=False)
    ax.set_title('XIST Expression by QC Status', fontsize=14, fontweight='bold')
    plt.tight_layout()
    return fig

def plot_highest_expr_genes(adata):
    """Top 20 expressed genes"""
    fig = plt.figure(figsize=(10, 6))
    sc.pl.highest_expr_genes(adata, n_top=20, show=False)
    return plt.gcf()

def plot_qc_metrics_density(adata):
    """QC metrics density plots"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_MT"]
    titles = ["Number of Genes", "Total Counts", "MT Percentage"]
    
    for ax, metric, title in zip(axes, metrics, titles):
        if metric in adata.obs.columns:
            sns.histplot(adata.obs[metric], ax=ax, kde=True, bins=100, stat='density')
            ax.set_xlabel(title, fontsize=12)
            ax.set_ylabel('Density', fontsize=12)
            ax.set_title(title, fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    return fig

def plot_counts_vs_genes(adata):
    """Scatter: total counts vs genes"""
    fig = plt.figure(figsize=(8, 6))
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_MT", show=False)
    return plt.gcf()

def plot_qclus_violins(adata):
    """Violin plots by QClus category"""
    desired_order = ['passed', 'scrublet filter', 'outlier filter', 'clustering filter', 'initial filter']
    qclus_order = [cat for cat in desired_order if cat in adata.obs['qclus'].unique()]
    
    adata.obs['qclus'] = adata.obs['qclus'].astype('category')
    adata.obs['qclus'] = adata.obs['qclus'].cat.reorder_categories(qclus_order, ordered=True)
    
    fig, axs = plt.subplots(2, 2, figsize=(15, 10))
    metrics = [
        ('fraction_unspliced', 'Fraction Unspliced'),
        ('pct_counts_MT', 'Pct Counts MT'),
        ('pct_counts_nuclear', 'Pct Counts Nuclear'),
        ('total_counts', 'Total Counts')
    ]
    
    for idx, (metric, title) in enumerate(metrics):
        row, col = idx // 2, idx % 2
        if metric in adata.obs.columns:
            sc.pl.violin(adata, metric, groupby='qclus', order=qclus_order, ax=axs[row, col], show=False)
            axs[row, col].set_title(title, fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    return fig

def plot_umap(adata_processed, title_prefix=""):
    """UMAP visualization"""
    metrics = ["fraction_unspliced", "pct_counts_MT", "pct_counts_nuclear", "qclus", 
               "pct_counts_CM_cyto", "pct_counts_CM_nucl", "pct_counts_FB", 
               "pct_counts_MP", "pct_counts_VEC", "pct_counts_SMC"]
    available_metrics = [m for m in metrics if m in adata_processed.obs.columns]
    sc.pl.umap(adata_processed, color=available_metrics, ncols=2, show=False)
    return plt.gcf()

def process_for_umap(adata):
    """Standard processing for UMAP"""
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.filter_genes(adata, min_cells=10)
    adata = adata[:, adata.var.highly_variable].copy()
    
    if 'pct_counts_MT' in adata.obs.columns:
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_MT'], n_jobs=4)
    else:
        sc.pp.regress_out(adata, ['total_counts'], n_jobs=4)
    
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='randomized')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added="leiden")
    return adata

def generate_html_report(sample_name, adata, output_dir):
    """Generate HTML report for one sample"""
    print(f"Processing {sample_name}...")
    
    # Basic stats
    n_cells = adata.n_obs
    n_genes = adata.n_vars
    qclus_counts = adata.obs['qclus'].value_counts().to_dict()
    qclus_table_html = adata.obs['qclus'].value_counts().to_frame(name='Count').to_html(classes='qclus-table')
    sex, _ = determine_sex(adata)
    
    # Generate plots
    print(f"  Generating plots...")
    img1 = fig_to_base64(plot_highest_expr_genes(adata))
    img2 = fig_to_base64(plot_qc_metrics_density(adata))
    img3 = fig_to_base64(plot_counts_vs_genes(adata))
    img4 = fig_to_base64(plot_qclus_violins(adata))
    img5 = fig_to_base64(plot_xist_violin(adata))
    
    # Process for UMAP
    print(f"  Processing for UMAP...")
    adata_viz = process_for_umap(adata.copy())
    img6 = fig_to_base64(plot_umap(adata_viz))
    
    # UMAP passed cells only
    adata_passed = adata_viz[adata_viz.obs.qclus == "passed"].copy()
    if adata_passed.n_obs > 0:
        sc.tl.pca(adata_passed, svd_solver='randomized')
        sc.pp.neighbors(adata_passed, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata_passed)
        sc.tl.leiden(adata_passed, key_added="leiden")
        img7 = fig_to_base64(plot_umap(adata_passed))
    else:
        fig7, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, 'No passed cells', ha='center', va='center', fontsize=14)
        ax.axis('off')
        img7 = fig_to_base64(fig7)
    
    # Generate HTML
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>QC Report - {sample_name}</title>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 20px; background-color: #f5f5f5; }}
        .container {{ max-width: 1400px; margin: 0 auto; background-color: white; padding: 30px; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 40px; border-left: 4px solid #3498db; padding-left: 15px; }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin: 30px 0; }}
        .stat-box {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 8px; text-align: center; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }}
        .stat-box.sex {{ background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); }}
        .stat-label {{ font-size: 14px; opacity: 0.9; margin-bottom: 8px; }}
        .stat-value {{ font-size: 32px; font-weight: bold; }}
        .plot-section {{ margin: 30px 0; padding: 20px; background-color: #fafafa; border-radius: 8px; }}
        .plot-section img {{ max-width: 100%; height: auto; display: block; margin: 20px auto; border-radius: 5px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }}
        .qclus-table {{ width: 100%; border-collapse: collapse; margin: 20px 0; font-size: 14px; }}
        .qclus-table th {{ background-color: #3498db; color: white; padding: 12px; text-align: left; }}
        .qclus-table td {{ padding: 10px; border-bottom: 1px solid #ddd; }}
        .qclus-table tr:hover {{ background-color: #f5f5f5; }}
        .timestamp {{ text-align: right; color: #7f8c8d; font-size: 12px; margin-top: 30px; }}
        .section-divider {{ height: 2px; background: linear-gradient(to right, #3498db, transparent); margin: 40px 0; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 snRNA-seq QC Report: {sample_name}</h1>
        <div class="stats-grid">
            <div class="stat-box"><div class="stat-label">Total Cells</div><div class="stat-value">{n_cells:,}</div></div>
            <div class="stat-box"><div class="stat-label">Total Genes</div><div class="stat-value">{n_genes:,}</div></div>
            <div class="stat-box"><div class="stat-label">Passed Cells</div><div class="stat-value">{qclus_counts.get('passed', 0):,}</div></div>
            <div class="stat-box sex"><div class="stat-label">Predicted Sex</div><div class="stat-value" style="font-size: 24px;">{sex}</div></div>
        </div>
        <h2>📊 QClus Filter Statistics</h2>
        <div class="plot-section">{qclus_table_html}</div>
        <div class="section-divider"></div>
        <h2>🧪 XIST Expression Analysis</h2>
        <div class="plot-section"><img src="{img5}" alt="XIST Expression"></div>
        <div class="section-divider"></div>
        <h2>📈 Highest Expressed Genes</h2>
        <div class="plot-section"><img src="{img1}" alt="Highest Expressed Genes"></div>
        <div class="section-divider"></div>
        <h2>📉 QC Metrics Distribution</h2>
        <div class="plot-section"><img src="{img2}" alt="QC Metrics Density"></div>
        <div class="section-divider"></div>
        <h2>🔬 Total Counts vs Genes</h2>
        <div class="plot-section"><img src="{img3}" alt="Counts vs Genes"></div>
        <div class="section-divider"></div>
        <h2>🎻 QC Metrics by Filter Status</h2>
        <div class="plot-section"><img src="{img4}" alt="QClus Violins"></div>
        <div class="section-divider"></div>
        <h2>🗺️ UMAP Visualization - All Cells</h2>
        <div class="plot-section"><img src="{img6}" alt="UMAP All Cells"></div>
        <div class="section-divider"></div>
        <h2>✅ UMAP Visualization - Passed Cells Only</h2>
        <div class="plot-section"><img src="{img7}" alt="UMAP Passed Cells"></div>
        <div class="timestamp">Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
    </div>
</body>
</html>"""
    
    output_file = os.path.join(output_dir, f"{sample_name}_QC_report.html")
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"  ✓ Report saved: {output_file}")
    
    return {
        'Sample': sample_name,
        'Total_Cells': n_cells,
        'Total_Genes': n_genes,
        'Passed_Cells': qclus_counts.get('passed', 0),
        'Pass_Rate': f"{qclus_counts.get('passed', 0)/n_cells*100:.1f}%",
        'Predicted_Sex': sex
    }

def main():
    if len(sys.argv) < 2:
        print("Usage: python qc_reports.py <base_directory>")
        print("Example: python qc_reports.py /mnt/archive/farhadie/public_studies/Bondareva_AF/cellranger_runs")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    
    # Find all h5ad files in subdirectories
    pattern = os.path.join(base_dir, "*/*/outs/*_QClus.h5ad")
    h5ad_files = glob.glob(pattern)
    
    if not h5ad_files:
        print(f"❌ No QClus h5ad files found in: {base_dir}")
        print(f"   Looking for pattern: */*/outs/*_QClus.h5ad")
        sys.exit(1)
    
    print(f"\n📁 Found {len(h5ad_files)} samples\n")
    
    # Create output directory
    output_dir = os.path.join(base_dir, "sample")
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}\n")
    
    # Process each sample
    summary_data = []
    for h5ad_file in h5ad_files:
        sample_name = os.path.basename(h5ad_file).replace('_QClus.h5ad', '')
        
        try:
            adata = sc.read_h5ad(h5ad_file)
            summary = generate_html_report(sample_name, adata, output_dir)
            summary_data.append(summary)
        except Exception as e:
            print(f"❌ Error processing {sample_name}: {str(e)}")
    
    # Create summary table
    summary_df = pd.DataFrame(summary_data)
    summary_csv = os.path.join(output_dir, "summary_table.csv")
    summary_df.to_csv(summary_csv, index=False)
    
    # Create summary HTML
    summary_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Summary - All Samples</title>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 20px; background-color: #f5f5f5; }}
        .container {{ max-width: 1400px; margin: 0 auto; background-color: white; padding: 30px; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        table {{ width: 100%; border-collapse: collapse; margin: 30px 0; font-size: 14px; }}
        th {{ background-color: #3498db; color: white; padding: 15px; text-align: left; font-weight: bold; }}
        td {{ padding: 12px 15px; border-bottom: 1px solid #ddd; }}
        tr:hover {{ background-color: #f5f5f5; }}
        tr:nth-child(even) {{ background-color: #fafafa; }}
        .sample-link {{ color: #3498db; text-decoration: none; font-weight: bold; }}
        .sample-link:hover {{ text-decoration: underline; }}
        .timestamp {{ text-align: right; color: #7f8c8d; font-size: 12px; margin-top: 30px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>📊 Summary Report - All Samples</h1>
        <table>
            <thead>
                <tr>
                    <th>Sample</th>
                    <th>Total Cells</th>
                    <th>Total Genes</th>
                    <th>Passed Cells</th>
                    <th>Pass Rate</th>
                    <th>Predicted Sex</th>
                    <th>Report</th>
                </tr>
            </thead>
            <tbody>"""
    
    for _, row in summary_df.iterrows():
        summary_html += f"""
                <tr>
                    <td>{row['Sample']}</td>
                    <td>{row['Total_Cells']}</td>
                    <td>{row['Total_Genes']}</td>
                    <td>{row['Passed_Cells']}</td>
                    <td>{row['Pass_Rate']}</td>
                    <td>{row['Predicted_Sex']}</td>
                    <td><a href="{row['Sample']}_QC_report.html" class="sample-link">View Report</a></td>
                </tr>"""
    
    summary_html += f"""
            </tbody>
        </table>
        <div class="timestamp">Summary generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
    </div>
</body>
</html>"""
    
    summary_html_file = os.path.join(output_dir, "summary.html")
    with open(summary_html_file, 'w', encoding='utf-8') as f:
        f.write(summary_html)
    
    print(f"\n✅ All reports generated!")
    print(f"📄 Summary CSV: {summary_csv}")
    print(f"🌐 Summary HTML: {summary_html_file}")
    print(f"\n📊 Summary:\n{summary_df.to_string(index=False)}")

if __name__ == "__main__":
    main()
