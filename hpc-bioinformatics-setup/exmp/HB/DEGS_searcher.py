import pandas as pd
import numpy as np
import json

# ============================================================================
# CONFIGURATION - Update these paths
# ============================================================================
DESEQ2_FEMALE = "DESeq2_results_CMF.csv"
DESEQ2_MALE = "DESeq2_results_CMM.csv"
OLGA_FEMALE = "Olga_CMF.csv"
OLGA_MALE = "Olga_CMM.csv"
OUTPUT_HTML = "gene_expression_report.html"

# Which layer to use for expression values
EXPRESSION_LAYER = 'normed_counts'  # or 'X', 'soupX_counts', 'psbulk_props'

# ============================================================================
# FUNCTIONS
# ============================================================================

def extract_gene_expression(gene_name, adata, use_layer='normed_counts'):
    """Extract gene expression values from AnnData object"""
    if gene_name not in adata.var_names:
        return None
    
    gene_idx = list(adata.var_names).index(gene_name)
    
    # Extract expression data
    if use_layer == 'X' or use_layer not in adata.layers:
        expr_data = adata.X[:, gene_idx]
    else:
        expr_data = adata.layers[use_layer][:, gene_idx]
    
    # Convert to array if sparse
    if hasattr(expr_data, 'toarray'):
        expr_data = expr_data.toarray().flatten()
    else:
        expr_data = np.array(expr_data).flatten()
    
    # Split by group
    control_mask = adata.obs['group'] == 'CT'
    patient_mask = adata.obs['group'] == 'AF'
    
    return {
        'control': expr_data[control_mask].tolist(),
        'patient': expr_data[patient_mask].tolist()
    }


def generate_report(cells_F, cells_M):
    """
    Generate HTML report from your cells_F and cells_M AnnData objects
    
    Parameters:
    -----------
    cells_F : AnnData
        Female cells (already filtered for cell type and sex)
    cells_M : AnnData
        Male cells (already filtered for cell type and sex)
    """
    
    print("=" * 70)
    print("Gene Expression Interactive Report Generator")
    print("=" * 70)
    
    # Load statistical results
    print("\nLoading statistical results...")
    deseq_f = pd.read_csv(DESEQ2_FEMALE, index_col=0)
    deseq_m = pd.read_csv(DESEQ2_MALE, index_col=0)
    olga_f = pd.read_csv(OLGA_FEMALE, sep=';', index_col=0)
    olga_m = pd.read_csv(OLGA_MALE, sep=';', index_col=0)
    
    print(f"✓ Loaded {len(deseq_f)} genes from DESeq2 Female")
    print(f"✓ Loaded {len(deseq_m)} genes from DESeq2 Male")
    print(f"✓ Loaded {len(olga_f)} genes from Olga Female")
    print(f"✓ Loaded {len(olga_m)} genes from Olga Male")
    
    # Get all unique genes
    all_genes = sorted(set(list(deseq_f.index) + list(deseq_m.index) + 
                           list(olga_f.index) + list(olga_m.index)))
    
    print(f"\nProcessing {len(all_genes)} unique genes...")
    
    # Prepare gene data
    gene_data = {}
    for i, gene in enumerate(all_genes, 1):
        if i % 1000 == 0:
            print(f"  Processed {i}/{len(all_genes)} genes...")
        
        # Get statistics
        stats = {}
        
        # DESeq2 Female
        if gene in deseq_f.index:
            padj = deseq_f.loc[gene, 'padj']
            lfc = deseq_f.loc[gene, 'log2FoldChange']
            stats['deseq_f_padj'] = float(padj) if pd.notna(padj) else None
            stats['deseq_f_lfc'] = float(lfc) if pd.notna(lfc) else None
        else:
            stats['deseq_f_padj'] = None
            stats['deseq_f_lfc'] = None
        
        # DESeq2 Male
        if gene in deseq_m.index:
            padj = deseq_m.loc[gene, 'padj']
            lfc = deseq_m.loc[gene, 'log2FoldChange']
            stats['deseq_m_padj'] = float(padj) if pd.notna(padj) else None
            stats['deseq_m_lfc'] = float(lfc) if pd.notna(lfc) else None
        else:
            stats['deseq_m_padj'] = None
            stats['deseq_m_lfc'] = None
        
        # Olga Female
        if gene in olga_f.index:
            padj = olga_f.loc[gene, 'CM_padj']
            lfc = olga_f.loc[gene, 'CM_log2FoldChange']
            stats['olga_f_padj'] = float(padj) if pd.notna(padj) else None
            stats['olga_f_lfc'] = float(lfc) if pd.notna(lfc) else None
        else:
            stats['olga_f_padj'] = None
            stats['olga_f_lfc'] = None
        
        # Olga Male
        if gene in olga_m.index:
            padj = olga_m.loc[gene, 'CM_padj']
            lfc = olga_m.loc[gene, 'CM_log2FoldChange']
            stats['olga_m_padj'] = float(padj) if pd.notna(padj) else None
            stats['olga_m_lfc'] = float(lfc) if pd.notna(lfc) else None
        else:
            stats['olga_m_padj'] = None
            stats['olga_m_lfc'] = None
        
        # Extract expression data
        female_expr = extract_gene_expression(gene, cells_F, EXPRESSION_LAYER)
        male_expr = extract_gene_expression(gene, cells_M, EXPRESSION_LAYER)
        
        # Store gene data
        gene_data[gene] = {
            'stats': stats,
            'boxplot': {
                'Female_Control': female_expr['control'] if female_expr else [],
                'Female_Patient': female_expr['patient'] if female_expr else [],
                'Male_Control': male_expr['control'] if male_expr else [],
                'Male_Patient': male_expr['patient'] if male_expr else []
            }
        }
    
    print(f"✓ Processed all {len(gene_data)} genes")
    
    # Generate HTML
    print("\nGenerating HTML report...")
    gene_data_json = json.dumps(gene_data)
    
    html_content = f'''<!DOCTYPE html>
<html>
<head>
    <title>Gene Expression Report</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        body {{
            font-family: Arial, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
        }}
        .search-box {{
            background: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        input[type="text"] {{
            width: 100%;
            padding: 15px;
            font-size: 16px;
            border: 2px solid #ddd;
            border-radius: 5px;
            box-sizing: border-box;
        }}
        input[type="text"]:focus {{
            outline: none;
            border-color: #667eea;
        }}
        .autocomplete {{
            position: relative;
        }}
        .autocomplete-items {{
            position: absolute;
            border: 1px solid #d4d4d4;
            border-bottom: none;
            border-top: none;
            z-index: 99;
            top: 100%;
            left: 0;
            right: 0;
            max-height: 300px;
            overflow-y: auto;
            background: white;
        }}
        .autocomplete-items div {{
            padding: 10px;
            cursor: pointer;
            background-color: #fff;
            border-bottom: 1px solid #d4d4d4;
        }}
        .autocomplete-items div:hover {{
            background-color: #e9e9e9;
        }}
        .results {{
            background: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            display: none;
        }}
        .results.show {{
            display: block;
        }}
        .gene-title {{
            font-size: 28px;
            color: #333;
            margin-bottom: 20px;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
        }}
        .plot-area {{
            margin: 30px 0;
            padding: 20px;
            background: #f9f9f9;
            border-radius: 5px;
        }}
        .stats {{
            margin-top: 30px;
            line-height: 2;
        }}
        .stats-line {{
            padding: 10px;
            background: #f9f9f9;
            margin: 5px 0;
            border-radius: 5px;
            font-size: 15px;
        }}
        .stats-label {{
            font-weight: bold;
            color: #555;
        }}
        .stats-value {{
            color: #667eea;
            font-weight: 600;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 Cardiocardiomyocytes Gene Expression Analysis </h1>
        <p>Pseudobulk performed by DecoupleR on 30 samples!!! </p>
    </div>
    
    <div class="search-box">
        <div class="autocomplete">
            <input type="text" id="geneInput" placeholder="Enter gene name (e.g., CDH11, IGFBP5, ADGRL3)..." />
            <div id="autocomplete-list" class="autocomplete-items"></div>
        </div>
    </div>
    
    <div id="results" class="results">
        <h2 class="gene-title" id="geneName"></h2>
        
        <div class="plot-area">
            <div id="boxplot"></div>
        </div>
        
        <div class="stats">
            <div class="stats-line">
                <span class="stats-label">The p adjusted value for female is:</span>
                <span class="stats-value" id="stat-deseq-f-padj"></span>
            </div>
            <div class="stats-line">
                <span class="stats-label">The log2 fold change for female is:</span>
                <span class="stats-value" id="stat-deseq-f-lfc"></span>
            </div>
            <div class="stats-line">
                <span class="stats-label">The p adjusted value for male is:</span>
                <span class="stats-value" id="stat-deseq-m-padj"></span>
            </div>
            <div class="stats-line">
                <span class="stats-label">The log2 fold change for male is:</span>
                <span class="stats-value" id="stat-deseq-m-lfc"></span>
            </div>
            <div class="stats-line">
                <span class="stats-label">The p adjusted value for female was reported based on Olga Analysis:</span>
                <span class="stats-value" id="stat-olga-f-padj"></span>
            </div>
            <div class="stats-line">
                <span class="stats-label">The log2 fold change for female was reported based on Olga Analysis:</span>
                <span class="stats-value" id="stat-olga-f-lfc"></span>
            </div>
            <div class="stats-line">
                <span class="stats-label">The p adjusted value for male was reported based on Olga Analysis:</span>
                <span class="stats-value" id="stat-olga-m-padj"></span>
            </div>
            <div class="stats-line">
                <span class="stats-label">The log2 fold change for male was reported based on Olga Analysis:</span>
                <span class="stats-value" id="stat-olga-m-lfc"></span>
            </div>
        </div>
    </div>

    <script>
        const geneData = {gene_data_json};
        const geneList = Object.keys(geneData).sort();
        
        const input = document.getElementById('geneInput');
        const autocompleteList = document.getElementById('autocomplete-list');
        const results = document.getElementById('results');
        
        // Autocomplete
        input.addEventListener('input', function() {{
            const val = this.value.trim().toUpperCase();
            closeAllLists();
            
            if (!val) return;
            
            const matches = geneList.filter(g => g.toUpperCase().includes(val)).slice(0, 10);
            
            matches.forEach(gene => {{
                const div = document.createElement('div');
                const idx = gene.toUpperCase().indexOf(val);
                div.innerHTML = gene.substr(0, idx) + 
                    '<strong>' + gene.substr(idx, val.length) + '</strong>' + 
                    gene.substr(idx + val.length);
                div.addEventListener('click', function() {{
                    input.value = gene;
                    closeAllLists();
                    showGene(gene);
                }});
                autocompleteList.appendChild(div);
            }});
        }});
        
        input.addEventListener('keydown', function(e) {{
            if (e.keyCode === 13) {{ // Enter key
                const gene = geneList.find(g => g.toUpperCase() === this.value.trim().toUpperCase());
                if (gene) showGene(gene);
            }}
        }});
        
        function closeAllLists() {{
            autocompleteList.innerHTML = '';
        }}
        
        document.addEventListener('click', function(e) {{
            if (e.target !== input) closeAllLists();
        }});
        
        function showGene(gene) {{
            if (!geneData[gene]) {{
                alert('Gene not found');
                return;
            }}
            
            const data = geneData[gene];
            
            // Update gene name
            document.getElementById('geneName').textContent = gene;
            
            // Update statistics
            document.getElementById('stat-deseq-f-padj').textContent = formatValue(data.stats.deseq_f_padj);
            document.getElementById('stat-deseq-f-lfc').textContent = formatValue(data.stats.deseq_f_lfc);
            document.getElementById('stat-deseq-m-padj').textContent = formatValue(data.stats.deseq_m_padj);
            document.getElementById('stat-deseq-m-lfc').textContent = formatValue(data.stats.deseq_m_lfc);
            document.getElementById('stat-olga-f-padj').textContent = formatValue(data.stats.olga_f_padj);
            document.getElementById('stat-olga-f-lfc').textContent = formatValue(data.stats.olga_f_lfc);
            document.getElementById('stat-olga-m-padj').textContent = formatValue(data.stats.olga_m_padj);
            document.getElementById('stat-olga-m-lfc').textContent = formatValue(data.stats.olga_m_lfc);
            
            // Create box plot
            const traces = [
                {{
                    y: data.boxplot.Female_Control,
                    name: 'Female Control',
                    type: 'box',
                    marker: {{ color: '#FFB6C1' }}
                }},
                {{
                    y: data.boxplot.Female_Patient,
                    name: 'Female Patient',
                    type: 'box',
                    marker: {{ color: '#DC143C' }}
                }},
                {{
                    y: data.boxplot.Male_Control,
                    name: 'Male Control',
                    type: 'box',
                    marker: {{ color: '#87CEEB' }}
                }},
                {{
                    y: data.boxplot.Male_Patient,
                    name: 'Male Patient',
                    type: 'box',
                    marker: {{ color: '#4169E1' }}
                }}
            ];
            
            const layout = {{
                title: `Expression Levels for ${{gene}}`,
                yaxis: {{ title: 'DESeq2 Normalized Counts' }},
                xaxis: {{ title: 'Group' }},
                height: 400
            }};
            
            Plotly.newPlot('boxplot', traces, layout);
            
            // Show results
            results.classList.add('show');
            results.scrollIntoView({{ behavior: 'smooth' }});
        }}
        
        function formatValue(val) {{
            if (val === null || val === undefined) return 'N/A';
            return val.toExponential(3);
        }}
    </script>
</body>
</html>'''
    
    # Save HTML file
    with open(OUTPUT_HTML, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"\n✓ Report successfully generated: {OUTPUT_HTML}")
    print(f"✓ Open the file in your web browser")
    print("=" * 70)


# ============================================================================
# USAGE
# ============================================================================

if __name__ == "__main__":
    # Check if cells_F and cells_M exist in the environment
    try:
        # Try to get from globals (when run in Jupyter/IPython)
        import __main__
        cells_F = getattr(__main__, 'cells_F', None)
        cells_M = getattr(__main__, 'cells_M', None)
        
        if cells_F is not None and cells_M is not None:
            print("Found cells_F and cells_M, generating report...")
            generate_report(cells_F, cells_M)
        else:
            print("=" * 70)
            print("Gene Expression Report Generator")
            print("=" * 70)
            print("\nERROR: cells_F and cells_M not found!")
            print("\nPlease run this script AFTER creating cells_F and cells_M:")
            print("\n  # Your code:")
            print("  cells_F = pdata[(pdata.obs['cell_type'] == 'CM') & (pdata.obs['sex'] == 'F')].copy()")
            print("  cells_M = pdata[(pdata.obs['cell_type'] == 'CM') & (pdata.obs['sex'] == 'M')].copy()")
            print("\n  # Then run:")
            print("  exec(open('gene_report.py').read())")
            print("\nOr use in Jupyter after defining cells_F and cells_M:")
            print("  %run gene_report.py")
            print("\n" + "=" * 70)
    except Exception as e:
        print(f"Error: {e}")
        print("\nTo use this script, ensure cells_F and cells_M are defined, then:")
        print("  from gene_report import generate_report")
        print("  generate_report(cells_F, cells_M)")
