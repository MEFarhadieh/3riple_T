#!/bin/bash
#SBATCH -J sNaM_modisco_1M
#SBATCH -p hpc
#SBATCH -t 122:00:00
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=256G
#SBATCH -o /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/logs/18_modisco_1M_%j.out
#SBATCH -e /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/logs/18_modisco_1M_%j.err

set -euo pipefail

export TMPDIR=/cl_tmp/farhadie/tmp
mkdir -p $TMPDIR

BASE=/mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted
SIF=/mnt/archive/farhadie/container/chrombpnet.sif
OUTDIR=$BASE/18_modisco_1M
mkdir -p $OUTDIR

echo "[$(date)] TF-MoDISco with 1M seqlets (CPU)"

if [ ! -f "$OUTDIR/ohe.npz" ]; then
    echo "ساخت ohe.npz و shap.npz..."
    singularity exec \
        --bind /mnt/archive --bind /cl_tmp \
        --env HDF5_PLUGIN_PATH=/opt/conda/lib/python3.9/site-packages/hdf5plugin/plugins \
        $SIF python3 << 'PYEOF'
import h5py
import numpy as np
BASE   = "/mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted"
OUTDIR = f"{BASE}/18_modisco_1M"
with h5py.File(f"{BASE}/06_contribution_scores/sNaM.profile_scores.h5", 'r') as f:
    ohe  = f['raw']['seq'][:].astype(np.float32)
    shap = f['shap']['seq'][:].astype(np.float32)
np.savez_compressed(f"{OUTDIR}/ohe.npz",  ohe)
np.savez_compressed(f"{OUTDIR}/shap.npz", shap)
print(f"saved: ohe={ohe.shape} shap={shap.shape}")
PYEOF
fi

echo "[$(date)] run modisco..."

singularity exec \
    --bind /mnt/archive --bind /cl_tmp \
    --env HDF5_PLUGIN_PATH=/opt/conda/lib/python3.9/site-packages/hdf5plugin/plugins \
    $SIF \
    modisco motifs \
        -s $OUTDIR/ohe.npz \
        -a $OUTDIR/shap.npz \
        -n 1000000 \
        -w 400 \
        -o $OUTDIR/modisco_results_1M.h5

echo "[$(date)] === Finished ==="
ls -lh $OUTDIR/

rm -f $OUTDIR/ohe.npz $OUTDIR/shap.npz
