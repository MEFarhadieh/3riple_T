# ================================
# HPC bioinformatics setup
# Add this to ~/.bashrc (Linux HPC)
# ================================

# Scratch and archive dirs (change if needed)
export LAB_SCRATCH="/cl_tmp/$USER"
export LAB_ARCHIVE="/mnt/archive/$USER"

# Nextflow work directory
export NXF_WORK="$LAB_SCRATCH/nf-work"
export NXF_HOME="$HOME/.nextflow"

# Temporary dirs (avoid filling $HOME)
export TMPDIR="$LAB_SCRATCH/tmp"
export TEMP="$LAB_SCRATCH/tmp"

# Singularity/Apptainer cache
export APPTAINER_CACHEDIR="$LAB_SCRATCH/singularity_cache"
export SINGULARITY_CACHEDIR="$LAB_SCRATCH/singularity_cache"
export NXF_SINGULARITY_CACHEDIR="$LAB_SCRATCH/singularity_cache"

# Conda/Mamba cache (optional)
export CONDA_PKGS_DIRS="$LAB_SCRATCH/conda_pkgs"
export MAMBA_ROOT_PREFIX="$LAB_SCRATCH/mamba"

# Auto-create directories if they don't exist
mkdir -p "$LAB_SCRATCH"/{nf-work,tmp,results,data,containers,logs} "$LAB_ARCHIVE"
