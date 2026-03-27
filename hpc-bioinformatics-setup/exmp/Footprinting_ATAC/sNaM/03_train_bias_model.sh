#!/bin/bash
#SBATCH -J sNaM_bias_train
#SBATCH -p gpu
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH -o /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/logs/03_bias_train_%j.out
#SBATCH -e /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/logs/03_bias_train_%j.err

set -euo pipefail

module load cuda/12.3

# ─────────────────────────────────────────────
# TMPDIR
# ─────────────────────────────────────────────
export TMPDIR=/cl_tmp/farhadie/tmp
export TEMP=/cl_tmp/farhadie/tmp
export TMP=/cl_tmp/farhadie/tmp
mkdir -p $TMPDIR

# ─────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────
BASE=/mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted
SIF=/mnt/archive/farhadie/container/chrombpnet.sif

BAM=$BASE/sNaM_sorted.bam
GENOME=/mnt/archive/farhadie/ref/mm10/mm10.fa
CHROMSIZES=$BASE/02_prep/mm10.main.chrom.sizes
PEAKS=$BASE/01_peaks/sNaM_peaks_no_blacklist_mainchrom.bed
NONPEAKS=$BASE/02_prep/sNaM_negatives.bed
FOLD=$BASE/02_prep/splits/fold_0.json

OUTDIR=$BASE/03_bias_model
mkdir -p $OUTDIR

echo "[$(date)] === Start training bias model ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# test GPU usage
singularity exec --nv --bind /mnt/archive --bind /cl_tmp $SIF \
    python3 -c "
import tensorflow as tf
gpus = tf.config.list_physical_devices('GPU')
print('GPUs found:', gpus)
if not gpus:
    raise RuntimeError('NO GPU DETECTED - aborting')
print('GPU OK, proceeding...')
"

singularity exec --nv \
    --bind /mnt/archive \
    --bind /cl_tmp \
    --env CUDA_VISIBLE_DEVICES=0 \
    --env TF_FORCE_GPU_ALLOW_GROWTH=true \
    $SIF \
    chrombpnet bias pipeline \
        -ibam $BAM \
        -d "ATAC" \
        -g $GENOME \
        -c $CHROMSIZES \
        -p $PEAKS \
        -n $NONPEAKS \
        -fl $FOLD \
        -b 0.5 \
        -o $OUTDIR

echo "[$(date)] === Bias model training is finished  ==="

if [ ! -f "$OUTDIR/models/bias.h5" ]; then
    echo "ERROR: file bias.h5 not found!" >&2
    exit 1
fi

echo " Model file path: $OUTDIR/models/bias.h5"
echo "Size: $(du -sh $OUTDIR/models/bias.h5)"
