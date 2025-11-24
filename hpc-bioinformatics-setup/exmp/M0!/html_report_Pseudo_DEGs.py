import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import numpy as np
import base64
from io import BytesIO
import os

# Define cell types
CELL_TYPES = ['AD', 'BL', 'CM', 'EEC', 'FB', 'MESO', 'MP', 'PER', 'SC', 'SMC', 'TL', 'VEC']
BASE_DIR = '/work/archive/public_studies/merged_data/DGE'


def fig_to_base64(fig):
    """Convert matplotlib figure to base64 string for embedding in HTML"""
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return img_base64


def create_venn_diagrams(ad_sig, female_sig, title_suffix):
    """Create Venn diagram and return as base64"""
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    venn = venn2([ad_sig, female_sig],
                 set_labels=('male', 'female'),
                 ax=ax)

    # Customize colors 
    if venn.get_patch_by_id('10'):
        venn.get_patch_by_id('10').set_color('#ADD8E6')  # Male
        venn.get_patch_by_id('10').set_alpha(0.8)
    if venn.get_patch_by_id('01'):
        venn.get_patch_by_id('01').set_color('#FFB6C1')  # Female
        venn.get_patch_by_id('01').set_alpha(0.8)
    if venn.get_patch_by_id('11'):
        venn.get_patch_by_id('11').set_color('#DDA0DD')  # Intersection
        venn.get_patch_by_id('11').set_alpha(0.9)

    ax.set_title(f'Significant Genes Intersection\n({title_suffix})',
                 fontsize=14, fontweight='bold', pad=20)

    return fig_to_base64(fig)


def create_volcano_plot(df, title):
    """
    Create a volcano plot (log2FC vs -log10(padj)) and return as base64 PNG.
    """
    required = {'log2FoldChange', 'padj'}
    if not required.issubset(df.columns):
        print(f"⚠ Volcano plot skipped for {title}: missing columns")
        return None

    tmp = df[['log2FoldChange', 'padj']].dropna()
    if tmp.empty:
        print(f"⚠ Volcano plot skipped for {title}: no non-NaN data")
        return None

    log2fc = tmp['log2FoldChange'].values
    padj = tmp['padj'].values
    neg_log10_padj = -np.log10(padj)

    sig_mask = padj < 0.05
    up_mask = sig_mask & (log2fc > 0)
    down_mask = sig_mask & (log2fc < 0)
    ns_mask = ~sig_mask

    fig, ax = plt.subplots(figsize=(7, 5))

    # Non-significant
    ax.scatter(
        log2fc[ns_mask], neg_log10_padj[ns_mask],
        s=6, alpha=0.4, edgecolors='none', color='#8B8989'
    )
    # Significant up
    ax.scatter(
        log2fc[up_mask], neg_log10_padj[up_mask],
        s=8, alpha=0.7, edgecolors='none', color='#FF0000'
    )
    # Significant down
    ax.scatter(
        log2fc[down_mask], neg_log10_padj[down_mask],
        s=8, alpha=0.7, edgecolors='none', color='#0000FF'
    )

    ax.axhline(-np.log10(0.05), color='grey', linestyle='--', linewidth=1)
    ax.axvline(0, color='grey', linestyle='--', linewidth=1)

    ax.set_xlabel('log2(Fold Change)')
    ax.set_ylabel('-log10(padj)')
    ax.set_title(f'Volcano Plot\n{title}')

    fig.tight_layout()
    return fig_to_base64(fig)


def analyze_cell_type(cell_type):
    """Analyze DESeq2 results for a given cell type"""

    # File paths
    male_file = os.path.join(BASE_DIR, cell_type, 'DESeq2_results_male.csv')
    female_file = os.path.join(BASE_DIR, cell_type, 'DESeq2_results_female.csv')

    # Check if files exist
    if not os.path.exists(male_file):
        print(f"⚠ Skipping {cell_type}: {male_file} not found")
        return None
    if not os.path.exists(female_file):
        print(f"⚠ Skipping {cell_type}: {female_file} not found")
        return None

    print(f"Processing {cell_type}...")

    # Read data
    try:
        male_results = pd.read_csv(male_file, index_col=0)
        female_results = pd.read_csv(female_file, index_col=0)
    except Exception as e:
        print(f"⚠ Error reading files for {cell_type}: {e}")
        return None

    # Get significant genes based on pvalue < 0.05
    male_sig_pvalue = set(male_results[male_results['pvalue'] < 0.05].index)
    female_sig_pvalue = set(female_results[female_results['pvalue'] < 0.05].index)

    # Get significant genes based on padj < 0.05
    male_sig_padj = set(male_results[male_results['padj'] < 0.05].index)
    female_sig_padj = set(female_results[female_results['padj'] < 0.05].index)

    # Calculate intersections
    intersection_pvalue = male_sig_pvalue & female_sig_pvalue
    intersection_padj = male_sig_padj & female_sig_padj

    # Create Venn diagrams
    venn_pvalue_img = create_venn_diagrams(male_sig_pvalue, female_sig_pvalue, 'p-value < 0.05')
    venn_padj_img = create_venn_diagrams(male_sig_padj, female_sig_padj, 'padj < 0.05')

    # Create volcano plots
    volcano_male_img = create_volcano_plot(male_results, f'{cell_type} - male')
    volcano_female_img = create_volcano_plot(female_results, f'{cell_type} - female')

    # Create intersection DataFrames
    intersection_df_pvalue = None
    intersection_df_padj = None

    if len(intersection_pvalue) > 0:
        intersection_df_pvalue = pd.DataFrame({
            'Gene': list(intersection_pvalue),
            'male_pvalue': [male_results.loc[g, 'pvalue'] for g in intersection_pvalue],
            'male_padj': [male_results.loc[g, 'padj'] for g in intersection_pvalue],
            'male_log2FC': [male_results.loc[g, 'log2FoldChange'] for g in intersection_pvalue],
            'male_baseMean': [male_results.loc[g, 'baseMean'] for g in intersection_pvalue],
            'female_pvalue': [female_results.loc[g, 'pvalue'] for g in intersection_pvalue],
            'female_padj': [female_results.loc[g, 'padj'] for g in intersection_pvalue],
            'female_log2FC': [female_results.loc[g, 'log2FoldChange'] for g in intersection_pvalue],
            'female_baseMean': [female_results.loc[g, 'baseMean'] for g in intersection_pvalue],
        }).sort_values('male_pvalue')

    if len(intersection_padj) > 0:
        intersection_df_padj = pd.DataFrame({
            'Gene': list(intersection_padj),
            'male_pvalue': [male_results.loc[g, 'pvalue'] for g in intersection_padj],
            'male_padj': [male_results.loc[g, 'padj'] for g in intersection_padj],
            'male_log2FC': [male_results.loc[g, 'log2FoldChange'] for g in intersection_padj],
            'male_baseMean': [male_results.loc[g, 'baseMean'] for g in intersection_padj],
            'female_pvalue': [female_results.loc[g, 'pvalue'] for g in intersection_padj],
            'female_padj': [female_results.loc[g, 'padj'] for g in intersection_padj],
            'female_log2FC': [female_results.loc[g, 'log2FoldChange'] for g in intersection_padj],
            'female_baseMean': [female_results.loc[g, 'baseMean'] for g in intersection_padj],
        }).sort_values('male_padj')

    # Create male-only significant genes DataFrames
    male_only_df_pvalue = None
    male_only_df_padj = None

    if len(male_sig_pvalue) > 0:
        male_only_df_pvalue = pd.DataFrame({
            'Gene': list(male_sig_pvalue),
            'pvalue': [male_results.loc[g, 'pvalue'] for g in male_sig_pvalue],
            'padj': [male_results.loc[g, 'padj'] for g in male_sig_pvalue],
            'log2FC': [male_results.loc[g, 'log2FoldChange'] for g in male_sig_pvalue],
            'baseMean': [male_results.loc[g, 'baseMean'] for g in male_sig_pvalue],
        }).sort_values('pvalue')

    if len(male_sig_padj) > 0:
        male_only_df_padj = pd.DataFrame({
            'Gene': list(male_sig_padj),
            'pvalue': [male_results.loc[g, 'pvalue'] for g in male_sig_padj],
            'padj': [male_results.loc[g, 'padj'] for g in male_sig_padj],
            'log2FC': [male_results.loc[g, 'log2FoldChange'] for g in male_sig_padj],
            'baseMean': [male_results.loc[g, 'baseMean'] for g in male_sig_padj],
        }).sort_values('padj')

    # Create female-only significant genes DataFrames
    female_only_df_pvalue = None
    female_only_df_padj = None

    if len(female_sig_pvalue) > 0:
        female_only_df_pvalue = pd.DataFrame({
            'Gene': list(female_sig_pvalue),
            'pvalue': [female_results.loc[g, 'pvalue'] for g in female_sig_pvalue],
            'padj': [female_results.loc[g, 'padj'] for g in female_sig_pvalue],
            'log2FC': [female_results.loc[g, 'log2FoldChange'] for g in female_sig_pvalue],
            'baseMean': [female_results.loc[g, 'baseMean'] for g in female_sig_pvalue],
        }).sort_values('pvalue')

    if len(female_sig_padj) > 0:
        female_only_df_padj = pd.DataFrame({
            'Gene': list(female_sig_padj),
            'pvalue': [female_results.loc[g, 'pvalue'] for g in female_sig_padj],
            'padj': [female_results.loc[g, 'padj'] for g in female_sig_padj],
            'log2FC': [female_results.loc[g, 'log2FoldChange'] for g in female_sig_padj],
            'baseMean': [female_results.loc[g, 'baseMean'] for g in female_sig_padj],
        }).sort_values('padj')

    return {
        'cell_type': cell_type,
        'male_shape': male_results.shape,
        'female_shape': female_results.shape,
        'male_sig_pvalue': len(male_sig_pvalue),
        'female_sig_pvalue': len(female_sig_pvalue),
        'intersection_pvalue': len(intersection_pvalue),
        'male_only_pvalue': len(male_sig_pvalue - female_sig_pvalue),
        'female_only_pvalue': len(female_sig_pvalue - male_sig_pvalue),
        'male_sig_padj': len(male_sig_padj),
        'female_sig_padj': len(female_sig_padj),
        'intersection_padj': len(intersection_padj),
        'male_only_padj': len(male_sig_padj - female_sig_padj),
        'female_only_padj': len(female_sig_padj - male_sig_padj),
        'venn_pvalue_img': venn_pvalue_img,
        'venn_padj_img': venn_padj_img,
        'intersection_df_pvalue': intersection_df_pvalue,
        'intersection_df_padj': intersection_df_padj,
        'male_df_pvalue': male_only_df_pvalue,
        'male_df_padj': male_only_df_padj,
        'female_df_pvalue': female_only_df_pvalue,
        'female_df_padj': female_only_df_padj,
        'volcano_male_img': volcano_male_img,
        'volcano_female_img': volcano_female_img,
    }


def dataframe_to_html_table(df, table_id):
    """Convert DataFrame to interactive HTML table with custom sorting and pagination"""
    if df is None or len(df) == 0:
        return '<p class="no-data">No significant genes in intersection</p>'

    import json

    print(f"  Creating table {table_id} with {len(df)} rows")

    # Convert DataFrame to list of dictionaries for JavaScript
    data_rows = []
    for idx in range(len(df)):
        row_data = {}
        for col in df.columns:
            if col == 'Gene':
                row_data[col] = str(df.iloc[idx][col])
            else:
                # Store both original value and formatted value
                original_val = float(df.iloc[idx][col])
                if 'pvalue' in col or 'padj' in col:
                    formatted_val = f'{original_val:.2e}'
                elif 'log2FC' in col:
                    formatted_val = f'{original_val:.3f}'
                elif 'baseMean' in col:
                    formatted_val = f'{original_val:.2f}'
                else:
                    formatted_val = str(original_val)
                row_data[col] = {'value': original_val, 'display': formatted_val}
        data_rows.append(row_data)

    # Convert to JSON string properly
    data_json = json.dumps(data_rows)
    print(f"  JSON length: {len(data_json)} characters")

    # Create HTML structure
    html = f'''
    <div class="table-controls" id="controls_{table_id}">
        <label>Show 
            <select onchange="changePageSize('{table_id}', this.value)">
                <option value="10" selected>10</option>
                <option value="50">50</option>
                <option value="all">All</option>
            </select>
        entries</label>
    </div>
    <div class="table-wrapper">
        <table id="{table_id}" class="interactive-table">
            <thead>
                <tr>
    '''

    # Add headers with click handlers
    for i, col in enumerate(df.columns):
        html += f'<th onclick="sortTable(\'{table_id}\', {i})" class="sortable">{col} <span class="sort-arrow">⇅</span></th>'

    html += f'''
                </tr>
            </thead>
            <tbody>
            </tbody>
        </table>
    </div>
    <div class="table-info" id="info_{table_id}"></div>
    <div class="pagination" id="pagination_{table_id}"></div>

    <script>
    window.tableStates['{table_id}'] = {{
        data: {data_json},
        sortColumn: -1,
        sortAscending: true,
        currentPage: 1,
        pageSize: 10
    }};
    renderTable('{table_id}');
    </script>
    '''

    return html


def generate_html_report(results):
    """Generate comprehensive HTML report"""

    html = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>DESeq2 Intersection Analysis - All Cell Types</title>
    
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            line-height: 1.6;
            color: #333;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
        }
        
        .container {
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 15px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow: hidden;
        }
        
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }
        
        .header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.2);
        }
        
        .header p {
            font-size: 1.2em;
            opacity: 0.9;
        }
        
        .nav {
            background: #f8f9fa;
            padding: 20px;
            border-bottom: 2px solid #e9ecef;
            position: sticky;
            top: 0;
            z-index: 100;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        
        .nav-buttons {
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            justify-content: center;
        }
        
        .nav-button {
            padding: 10px 20px;
            background: #667eea;
            color: white;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            font-size: 0.95em;
            font-weight: 600;
            transition: all 0.3s ease;
            text-decoration: none;
        }
        
        .nav-button:hover {
            background: #764ba2;
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(102, 126, 234, 0.4);
        }
        
        .summary {
            padding: 40px;
            background: #f8f9fa;
        }
        
        .summary h2 {
            color: #667eea;
            margin-bottom: 20px;
            font-size: 2em;
        }
        
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }
        
        .summary-card {
            background: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            text-align: center;
            transition: transform 0.3s ease;
        }
        
        .summary-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 5px 20px rgba(0,0,0,0.15);
        }
        
        .summary-card h3 {
            color: #667eea;
            font-size: 1.2em;
            margin-bottom: 10px;
        }
        
        .summary-card .number {
            font-size: 2.5em;
            font-weight: bold;
            color: #764ba2;
        }
        
        .cell-type-section {
            padding: 40px;
            border-bottom: 3px solid #e9ecef;
        }
        
        .cell-type-section:nth-child(even) {
            background: #f8f9fa;
        }
        
        .cell-type-header {
            display: flex;
            align-items: center;
            margin-bottom: 30px;
            padding-bottom: 15px;
            border-bottom: 3px solid #667eea;
        }
        
        .cell-type-header h2 {
            color: #667eea;
            font-size: 2.2em;
            margin-right: 15px;
        }
        
        .cell-type-badge {
            background: #667eea;
            color: white;
            padding: 5px 15px;
            border-radius: 20px;
            font-size: 0.9em;
            font-weight: 600;
        }
        
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin: 30px 0;
        }
        
        .stat-box {
            background: white;
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        
        .stat-box h4 {
            color: #666;
            font-size: 0.9em;
            margin-bottom: 5px;
        }
        
        .stat-box .value {
            font-size: 1.8em;
            font-weight: bold;
            color: #667eea;
        }
        
        .venn-container {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 30px;
            margin: 30px 0;
        }
        
        .venn-box {
            background: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            text-align: center;
        }
        
        .venn-box h3 {
            color: #667eea;
            margin-bottom: 15px;
            font-size: 1.3em;
        }
        
        .venn-box img {
            max-width: 100%;
            height: auto;
            border-radius: 8px;
        }

        /* Volcano plot layout */
        .volcano-container {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 30px;
            margin: 30px 0;
        }
        
        .volcano-box {
            background: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            text-align: center;
        }
        
        .volcano-box h3 {
            color: #667eea;
            margin-bottom: 15px;
            font-size: 1.3em;
        }
        
        .volcano-box img {
            max-width: 100%;
            height: auto;
            border-radius: 8px;
        }
        
        .table-section {
            margin-top: 40px;
        }
        
        .table-section h3 {
            color: #667eea;
            margin-bottom: 20px;
            font-size: 1.5em;
            display: flex;
            align-items: center;
        }
        
        .table-section h3::before {
            content: "📊";
            margin-right: 10px;
        }
        
        .table-controls {
            background: white;
            padding: 15px 20px;
            border-radius: 10px 10px 0 0;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            display: flex;
            justify-content: flex-start;
            align-items: center;
            margin-bottom: 0;
        }
        
        .table-controls label {
            font-size: 0.95em;
            color: #666;
            display: flex;
            align-items: center;
            gap: 10px;
        }
        
        .table-controls select {
            padding: 8px 12px;
            border: 2px solid #667eea;
            border-radius: 6px;
            background: white;
            color: #667eea;
            font-size: 0.95em;
            font-weight: 600;
            cursor: pointer;
            outline: none;
            transition: all 0.3s ease;
        }
        
        .table-controls select:hover {
            background: #f0f4ff;
            border-color: #764ba2;
        }
        
        .table-controls select:focus {
            border-color: #764ba2;
            box-shadow: 0 0 0 3px rgba(118, 75, 162, 0.1);
        }
        
        .table-wrapper {
            background: white;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            overflow-x: auto;
            margin-bottom: 0;
        }
        
        .interactive-table {
            width: 100%;
            border-collapse: collapse;
            font-size: 0.9em;
        }
        
        .interactive-table thead th {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
            cursor: pointer;
            user-select: none;
            position: relative;
            border-bottom: 3px solid #764ba2;
        }
        
        .interactive-table thead th:hover {
            background: linear-gradient(135deg, #5568d3 0%, #654a8e 100%);
        }
        
        .interactive-table thead th .sort-arrow {
            margin-left: 5px;
            opacity: 0.5;
            font-size: 0.9em;
        }
        
        .interactive-table thead th.sorted .sort-arrow {
            opacity: 1;
        }
        
        .interactive-table tbody td {
            padding: 12px 15px;
            border-bottom: 1px solid #e9ecef;
        }
        
        .interactive-table tbody tr:hover {
            background: #f0f4ff;
        }
        
        .interactive-table tbody tr:nth-child(even) {
            background: #fafbfc;
        }
        
        .interactive-table tbody tr:nth-child(even):hover {
            background: #f0f4ff;
        }
        
        .table-info {
            background: white;
            padding: 12px 20px;
            border-top: 1px solid #e9ecef;
            font-size: 0.9em;
            color: #666;
        }
        
        .pagination {
            background: white;
            padding: 15px 20px;
            border-radius: 0 0 10px 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            display: flex;
            justify-content: center;
            align-items: center;
            gap: 10px;
            margin-bottom: 30px;
        }
        
        .pagination button {
            padding: 8px 16px;
            border: 2px solid #667eea;
            border-radius: 6px;
            background: white;
            color: #667eea;
            font-size: 0.9em;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.3s ease;
        }
        
        .pagination button:hover:not(:disabled) {
            background: #667eea;
            color: white;
        }
        
        .pagination button:disabled {
            opacity: 0.3;
            cursor: not-allowed;
        }
        
        .pagination .page-numbers {
            display: flex;
            gap: 5px;
        }
        
        .pagination .page-number {
            padding: 8px 12px;
            border: 2px solid #667eea;
            border-radius: 6px;
            background: white;
            color: #667eea;
            font-size: 0.9em;
            font-weight: 600;
            cursor: pointer;
            min-width: 40px;
            text-align: center;
            transition: all 0.3s ease;
        }
        
        .pagination .page-number:hover {
            background: #f0f4ff;
        }
        
        .pagination .page-number.active {
            background: #667eea;
            color: white;
        }
        
        .no-data {
            text-align: center;
            padding: 40px;
            color: #999;
            font-style: italic;
        }
        
        .footer {
            background: #2d3748;
            color: white;
            padding: 30px;
            text-align: center;
        }
        
        .footer p {
            margin: 5px 0;
            opacity: 0.8;
        }
        
        @media (max-width: 768px) {
            .header h1 {
                font-size: 1.8em;
            }
            
            .venn-container, .volcano-container {
                grid-template-columns: 1fr;
            }
            
            .nav-buttons {
                flex-direction: column;
            }
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 DESeq2 Pseudo-bulk Analysis</h1>
            <p> Analysis Across All Cell Types</p>
        </div>
        
        <div class="nav">
            <div class="nav-buttons">
                <a href="#summary" class="nav-button">📊 Summary</a>
"""

    # Add navigation buttons for each cell type
    for result in results:
        if result is not None:
            html += f'                <a href="#{result["cell_type"]}" class="nav-button">{result["cell_type"]}</a>\n'

    html += """
            </div>
        </div>
        
        <script>
        // DEFINE ALL FUNCTIONS HERE FIRST - BEFORE ANY TABLES
        window.tableStates = {};
        
        function renderTable(tableId) {
            const state = window.tableStates[tableId];
            if (!state) return;
            
            const table = document.getElementById(tableId);
            if (!table) return;
            
            const tbody = table.querySelector('tbody');
            
            let sortedData = [...state.data];
            if (state.sortColumn >= 0) {
                const colName = Object.keys(state.data[0])[state.sortColumn];
                sortedData.sort((a, b) => {
                    let aVal, bVal;
                    if (colName === 'Gene') {
                        aVal = a[colName];
                        bVal = b[colName];
                        return state.sortAscending ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
                    } else {
                        aVal = a[colName].value;
                        bVal = b[colName].value;
                        return state.sortAscending ? aVal - bVal : bVal - aVal;
                    }
                });
            }
            
            const totalRows = sortedData.length;
            const pageSize = state.pageSize === 'all' ? totalRows : parseInt(state.pageSize);
            const startIdx = (state.currentPage - 1) * pageSize;
            const endIdx = Math.min(startIdx + pageSize, totalRows);
            const pageData = sortedData.slice(startIdx, endIdx);
            
            tbody.innerHTML = '';
            pageData.forEach(row => {
                const tr = document.createElement('tr');
                Object.keys(row).forEach(key => {
                    const td = document.createElement('td');
                    td.textContent = key === 'Gene' ? row[key] : row[key].display;
                    tr.appendChild(td);
                });
                tbody.appendChild(tr);
            });
            
            const headers = table.querySelectorAll('thead th');
            headers.forEach((th, idx) => {
                th.classList.remove('sorted');
                const arrow = th.querySelector('.sort-arrow');
                if (arrow) {
                    if (idx === state.sortColumn) {
                        th.classList.add('sorted');
                        arrow.textContent = state.sortAscending ? '↑' : '↓';
                    } else {
                        arrow.textContent = '⇅';
                    }
                }
            });
            
            const infoDiv = document.getElementById('info_' + tableId);
            if (infoDiv) {
                infoDiv.textContent = 'Showing ' + (startIdx + 1) + ' to ' + endIdx + ' of ' + totalRows + ' genes';
            }
            
            renderPagination(tableId, state.currentPage, Math.ceil(totalRows / pageSize));
        }
        
        function renderPagination(tableId, currentPage, totalPages) {
            const paginationDiv = document.getElementById('pagination_' + tableId);
            if (!paginationDiv || totalPages <= 1) {
                if (paginationDiv) paginationDiv.innerHTML = '';
                return;
            }
            
            let html = '<button onclick="changePage(\\'' + tableId + '\\', ' + (currentPage - 1) + ')" ' + 
                       (currentPage === 1 ? 'disabled' : '') + '>Previous</button>';
            html += '<div class="page-numbers">';
            
            for (let i = 1; i <= Math.min(totalPages, 10); i++) {
                html += '<div class="page-number ' + (i === currentPage ? 'active' : '') + '" ' +
                        'onclick="changePage(\\'' + tableId + '\\', ' + i + ')">' + i + '</div>';
            }
            
            html += '</div>';
            html += '<button onclick="changePage(\\'' + tableId + '\\', ' + (currentPage + 1) + ')" ' +
                    (currentPage === totalPages ? 'disabled' : '') + '>Next</button>';
            
            paginationDiv.innerHTML = html;
        }
        
        function sortTable(tableId, columnIndex) {
            const state = window.tableStates[tableId];
            if (!state) return;
            if (state.sortColumn === columnIndex) {
                state.sortAscending = !state.sortAscending;
            } else {
                state.sortColumn = columnIndex;
                state.sortAscending = true;
            }
            renderTable(tableId);
        }
        
        function changePage(tableId, newPage) {
            const state = window.tableStates[tableId];
            if (!state) return;
            state.currentPage = newPage;
            renderTable(tableId);
        }
        
        function changePageSize(tableId, newSize) {
            const state = window.tableStates[tableId];
            if (!state) return;
            state.pageSize = newSize;
            state.currentPage = 1;
            renderTable(tableId);
        }
        </script>
        
        <div id="summary" class="summary">
            <h2>📊 Analysis Summary</h2>
            <div class="summary-grid">
"""

    # Add summary statistics
    total_cell_types = len([r for r in results if r is not None])
    total_intersection_pvalue = sum([r['intersection_pvalue'] for r in results if r is not None])
    total_intersection_padj = sum([r['intersection_padj'] for r in results if r is not None])

    html += f"""
                <div class="summary-card">
                    <h3>Cell Types Analyzed</h3>
                    <div class="number">{total_cell_types}</div>
                </div>
                <div class="summary-card">
                    <h3>Total Genes (p-value)</h3>
                    <div class="number">{total_intersection_pvalue}</div>
                </div>
                <div class="summary-card">
                    <h3>Total Genes (padj)</h3>
                    <div class="number">{total_intersection_padj}</div>
                </div>
"""

    html += """
            </div>
        </div>
"""

    # Add sections for each cell type
    for result in results:
        if result is None:
            continue

        ct = result['cell_type']

        html += f"""
        <div id="{ct}" class="cell-type-section">
            <div class="cell-type-header">
                <h2>{ct}</h2>
                <span class="cell-type-badge">Cell Type</span>
            </div>
            
            <div class="stats-grid">
                <div class="stat-box">
                    <h4>male Total Genes</h4>
                    <div class="value">{result['male_shape'][0]}</div>
                </div>
                <div class="stat-box">
                    <h4>female Total Genes</h4>
                    <div class="value">{result['female_shape'][0]}</div>
                </div>
                <div class="stat-box">
                    <h4>male Sig (p-value)</h4>
                    <div class="value">{result['male_sig_pvalue']}</div>
                </div>
                <div class="stat-box">
                    <h4>female Sig (p-value)</h4>
                    <div class="value">{result['female_sig_pvalue']}</div>
                </div>
                <div class="stat-box">
                    <h4>Intersection (p-value)</h4>
                    <div class="value">{result['intersection_pvalue']}</div>
                </div>
                <div class="stat-box">
                    <h4>male Sig (padj)</h4>
                    <div class="value">{result['male_sig_padj']}</div>
                </div>
                <div class="stat-box">
                    <h4>female Sig (padj)</h4>
                    <div class="value">{result['female_sig_padj']}</div>
                </div>
                <div class="stat-box">
                    <h4>Intersection (padj)</h4>
                    <div class="value">{result['intersection_padj']}</div>
                </div>
            </div>
            
            <div class="venn-container">
                <div class="venn-box">
                    <h3>P-value &lt; 0.05</h3>
                    <img src="data:image/png;base64,{result['venn_pvalue_img']}" alt="Venn diagram p-value">
                </div>
                <div class="venn-box">
                    <h3>Adjusted P-value &lt; 0.05</h3>
                    <img src="data:image/png;base64,{result['venn_padj_img']}" alt="Venn diagram padj">
                </div>
            </div>

            <div class="volcano-container">
                <div class="volcano-box">
                    <h3>Volcano Plot – male</h3>
                    {(
                        f'<img src="data:image/png;base64,{result["volcano_male_img"]}" alt="Volcano male">'
                        if result["volcano_male_img"]
                        else '<p class="no-data">Volcano plot not available</p>'
                    )}
                </div>
                <div class="volcano-box">
                    <h3>Volcano Plot – female</h3>
                    {(
                        f'<img src="data:image/png;base64,{result["volcano_female_img"]}" alt="Volcano female">'
                        if result["volcano_female_img"]
                        else '<p class="no-data">Volcano plot not available</p>'
                    )}
                </div>
            </div>
            
            <div class="table-section">
                <h3>Male Significant Genes (p-value &lt; 0.05)</h3>
                {dataframe_to_html_table(result['male_df_pvalue'], f'table_{ct}_male_pvalue')}
            </div>
            
            <div class="table-section">
                <h3>Male Significant Genes (padj &lt; 0.05)</h3>
                {dataframe_to_html_table(result['male_df_padj'], f'table_{ct}_male_padj')}
            </div>
            
            <div class="table-section">
                <h3>Female Significant Genes (p-value &lt; 0.05)</h3>
                {dataframe_to_html_table(result['female_df_pvalue'], f'table_{ct}_female_pvalue')}
            </div>
            
            <div class="table-section">
                <h3>Female Significant Genes (padj &lt; 0.05)</h3>
                {dataframe_to_html_table(result['female_df_padj'], f'table_{ct}_female_padj')}
            </div>
            
            <div class="table-section">
                <h3>Intersection - Significant Genes (p-value &lt; 0.05)</h3>
                {dataframe_to_html_table(result['intersection_df_pvalue'], f'table_{ct}_pvalue')}
            </div>
            
            <div class="table-section">
                <h3>Intersection - Significant Genes (padj &lt; 0.05)</h3>
                {dataframe_to_html_table(result['intersection_df_padj'], f'table_{ct}_padj')}
            </div>
        </div>
"""

    html += """
        <div class="footer">
            <p><strong>DESeq2 Intersection Analysis Report</strong></p>
            <p>Generated with Python • matplotlib-venn • pandas</p>
            <p>male = DESeq2_results_male.csv | female = DESeq2_results_female.csv</p>
        </div>
    </div>
    
    <script>
        // Smooth scrolling
        document.querySelectorAll('a[href^="#"]').forEach(anchor => {
            anchor.addEventListener('click', function (e) {
                e.preventDefault();
                const target = document.querySelector(this.getAttribute('href'));
                if (target) {
                    target.scrollIntoView({ behavior: 'smooth', block: 'start' });
                }
            });
        });
    </script>
</body>
</html>
"""

    return html


# Main execution
print("=" * 70)
print("DESeq2 INTERSECTION ANALYSIS - ALL CELL TYPES")
print("=" * 70)

results = []
for cell_type in CELL_TYPES:
    result = analyze_cell_type(cell_type)
    results.append(result)
    if result:
        print(f"✓ {cell_type}: {result['intersection_pvalue']} genes (p-value), {result['intersection_padj']} genes (padj)")

print("\n" + "=" * 70)
print("GENERATING HTML REPORT")
print("=" * 70)

html_content = generate_html_report(results)

output_file = 'deseq2_female_vs_male_report.html'
with open(output_file, 'w', encoding='utf-8') as f:
    f.write(html_content)

print(f"✓ HTML report generated: {output_file}")
print("=" * 70)
print("ANALYSIS COMPLETE!")
print("=" * 70)
print(f"\nProcessed {len([r for r in results if r is not None])} cell types")
print(f"Output: {output_file}")
