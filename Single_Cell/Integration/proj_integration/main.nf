nextflow.enable.dsl=2

params.input     = "meta/samplesheet.csv"
params.outdir    = "results"
params.batch_key = "batch"
params.sex_key   = "sex"
params.methods   = "seurat,harmony,scvi,liger"   
params.hvgs      = "2000,3000"
params.pcs       = "30,50"
params.neighbors = "30"

process PREP {
  conda "envs/py_scvi.yml"
  input:
    path samplesheet from file(params.input)
  output:
    path "data/combined_raw.h5ad"
  script:
  """
  python bin/prep_inputs.py
  """
}

def grid = Channel
  .from(params.hvgs.split(','))
  .combine(Channel.from(params.pcs.split(',')))
  .combine(Channel.from(params.neighbors.split(',')))

process HARMONY {
  tag "${hvgs}_${pcs}_${nn}"
  conda "envs/py_scvi.yml"
  input:
    path adata from PREP.out
    val hvgs from grid.map{ it[0] }
    val pcs  from grid.map{ it[1] }
    val nn   from grid.map{ it[2] }
  when:
    params.methods.contains("harmony")
  output:
    path "results/run_h${hvgs}_p${pcs}_k${nn}/harmony/integrated.h5ad"
    path "results/run_h${hvgs}_p${pcs}_k${nn}/harmony/umap.png"
  script:
  """
  mkdir -p results/run_h${hvgs}_p${pcs}_k${nn}/harmony
  python bin/integrate_harmony.py \
    --infile data/combined_raw.h5ad \
    --outfile results/run_h${hvgs}_p${pcs}_k${nn}/harmony/integrated.h5ad \
    --batch_key ${params.batch_key} --n_hvgs ${hvgs} --n_pcs ${pcs} --neighbors ${nn}
  python bin/plot_umap_panels.py \
    --infile results/run_h${hvgs}_p${pcs}_k${nn}/harmony/integrated.h5ad \
    --outfile results/run_h${hvgs}_p${pcs}_k${nn}/harmony/umap.png
  """
}

process SCVI {
  tag "${hvgs}_${pcs}_${nn}"
  conda "envs/py_scvi.yml"
  input:
    path adata from PREP.out
    val hvgs from grid.map{ it[0] }
    val pcs  from grid.map{ it[1] }  
    val nn   from grid.map{ it[2] }
  when:
    params.methods.contains("scvi")
  output:
    path "results/run_h${hvgs}_p${pcs}_k${nn}/scvi/integrated.h5ad"
    path "results/run_h${hvgs}_p${pcs}_k${nn}/scvi/umap.png"
  script:
  """
  mkdir -p results/run_h${hvgs}_p${pcs}_k${nn}/scvi
  python bin/integrate_scvi.py \
    --infile data/combined_raw.h5ad \
    --outfile results/run_h${hvgs}_p${pcs}_k${nn}/scvi/integrated.h5ad \
    --batch_key ${params.batch_key} --n_hvgs ${hvgs} --neighbors ${nn}
  python bin/plot_umap_panels.py \
    --infile results/run_h${hvgs}_p${pcs}_k${nn}/scvi/integrated.h5ad \
    --outfile results/run_h${hvgs}_p${pcs}_k${nn}/scvi/umap.png
  """
}

process SEURAT_RPCA {
  tag "${hvgs}_${pcs}_${nn}"
  conda "envs/r_seurat.yml"
  input:
    path adata from PREP.out
    val hvgs from grid.map{ it[0] }
    val pcs  from grid.map{ it[1] }
    val nn   from grid.map{ it[2] }
  when:
    params.methods.contains("seurat")
  output:
    path "results/run_h${hvgs}_p${pcs}_k${nn}/seurat/integrated.h5ad"
    path "results/run_h${hvgs}_p${pcs}_k${nn}/seurat/umap.png"
  script:
  """
  mkdir -p results/run_h${hvgs}_p${pcs}_k${nn}/seurat
  Rscript bin/integrate_seurat_rpca.R \
    --infile data/combined_raw.h5ad \
    --outfile results/run_h${hvgs}_p${pcs}_k${nn}/seurat/integrated.h5ad \
    --batch_key ${params.batch_key} --n_hvgs ${hvgs} --n_pcs ${pcs} --neighbors ${nn}
  python bin/plot_umap_panels.py \
    --infile results/run_h${hvgs}_p${pcs}_k${nn}/seurat/integrated.h5ad \
    --outfile results/run_h${hvgs}_p${pcs}_k${nn}/seurat/umap.png
  """
}

process LIGER {
  tag "${hvgs}_${pcs}_${nn}"
  conda "envs/r_seurat.yml"
  input:
    path adata from PREP.out
    val hvgs from grid.map{ it[0] }
    val pcs  from grid.map{ it[1] }
    val nn   from grid.map{ it[2] }
  when:
    params.methods.contains("liger")
  output:
    path "results/run_h${hvgs}_p${pcs}_k${nn}/liger/integrated.h5ad"
    path "results/run_h${hvgs}_p${pcs}_k${nn}/liger/umap.png"
  script:
  """
  mkdir -p results/run_h${hvgs}_p${pcs}_k${nn}/liger
  Rscript bin/integrate_liger.R \
    --infile data/combined_raw.h5ad \
    --outfile results/run_h${hvgs}_p${pcs}_k${nn}/liger/integrated.h5ad \
    --batch_key ${params.batch_key} --neighbors ${nn}
  python bin/plot_umap_panels.py \
    --infile results/run_h${hvgs}_p${pcs}_k${nn}/liger/integrated.h5ad \
    --outfile results/run_h${hvgs}_p${pcs}_k${nn}/liger/umap.png
  """
}

process METRICS {
  tag "${method}_${hvgs}_${pcs}_${nn}"
  conda "envs/py_scvi.yml"
  input:
    path infile
    val method
    val hvgs
    val pcs
    val nn
  output:
    path "results/run_h${hvgs}_p${pcs}_k${nn}/${method}/metrics.csv"
  script:
  def embed_map = [
    "harmony":"X_harmony",
    "scvi":"X_scvi",
    "seurat":"X_pca",     
    "liger":"X_liger"     
  ]
  """
  python bin/metrics_scib.py \
    --infile ${infile} \
    --embed_key ${embed_map[method]} \
    --batch_key ${params.batch_key} \
    --sex_key ${params.sex_key} \
    --outfile results/run_h${hvgs}_p${pcs}_k${nn}/${method}/metrics.csv
  """
}

workflow {
  Channel.fromPath("results/run_h*_p*_k*/*/integrated.h5ad").ifEmpty([]).set { PREV }  
  def runs = []

  if (params.methods.contains("harmony")) {
    HARMONY(PREP.out, grid.map{ it[0] }, grid.map{ it[1] }, grid.map{ it[2] })
      .view()
    runs << Channel.fromPath("results/run_h*_p*_k*/harmony/integrated.h5ad").map{ [it,"harmony", it.parent.name.split("_")[1].substring(1), it.parent.name.split("_")[2].substring(1), it.parent.name.split("_")[3].substring(1)] }
  }
  if (params.methods.contains("scvi")) {
    SCVI(PREP.out, grid.map{ it[0] }, grid.map{ it[1] }, grid.map{ it[2] })
    runs << Channel.fromPath("results/run_h*_p*_k*/scvi/integrated.h5ad").map{ [it,"scvi", it.parent.name.split("_")[1].substring(1), it.parent.name.split("_")[2].substring(1), it.parent.name.split("_")[3].substring(1)] }
  }
  if (params.methods.contains("seurat")) {
    SEURAT_RPCA(PREP.out, grid.map{ it[0] }, grid.map{ it[1] }, grid.map{ it[2] })
    runs << Channel.fromPath("results/run_h*_p*_k*/seurat/integrated.h5ad").map{ [it,"seurat", it.parent.name.split("_")[1].substring(1), it.parent.name.split("_")[2].substring(1), it.parent.name.split("_")[3].substring(1)] }
  }
  if (params.methods.contains("liger")) {
    LIGER(PREP.out, grid.map{ it[0] }, grid.map{ it[1] }, grid.map{ it[2] })
    runs << Channel.fromPath("results/run_h*_p*_k*/liger/integrated.h5ad").map{ [it,"liger", it.parent.name.split("_")[1].substring(1), it.parent.name.split("_")[2].substring(1), it.parent.name.split("_")[3].substring(1)] }
  }

  merged = Channel.concat(runs).flatten().set { TO_METRICS }
  METRICS( TO_METRICS.map{ it[0] }, TO_METRICS.map{ it[1] }, TO_METRICS.map{ it[2] }, TO_METRICS.map{ it[3] }, TO_METRICS.map{ it[4] } )
}
