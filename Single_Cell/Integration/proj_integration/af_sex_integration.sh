################################################################################################
# 0) Build the directories structure
################################################################################################
cd "$LAB_SCRATCH/projects"
mkdir -p $LAB_SCRATCH/projects/af_sex_integration/{data,meta,envs,bin,report}

touch $LAB_SCRATCH/projects/af_sex_integration/meta/samplesheet.csv
touch $LAB_SCRATCH/projects/af_sex_integration/envs/py_scvi.yml
touch $LAB_SCRATCH/projects/af_sex_integration/envs/r_seurat.yml
touch $LAB_SCRATCH/projects/af_sex_integration/bin/prep_inputs.py
touch $LAB_SCRATCH/projects/af_sex_integration/bin/integrate_harmony.py
touch $LAB_SCRATCH/projects/af_sex_integration/bin/integrate_scvi.py
touch $LAB_SCRATCH/projects/af_sex_integration/bin/integrate_seurat_rpca.R
touch $LAB_SCRATCH/projects/af_sex_integration/bin/integrate_liger.R
touch $LAB_SCRATCH/projects/af_sex_integration/bin/metrics_scib.py
touch $LAB_SCRATCH/projects/af_sex_integration/bin/plot_umap_panels.py
touch $LAB_SCRATCH/projects/af_sex_integration/report/report.qmd
touch $LAB_SCRATCH/projects/af_sex_integration/main.nf
touch $LAB_SCRATCH/projects/af_sex_integration/nextflow.config

tree -lh
