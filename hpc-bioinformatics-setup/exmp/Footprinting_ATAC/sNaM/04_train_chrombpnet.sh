#!/bin/bash
#SBATCH -J sNaM_chrombpnet_train
#SBATCH -p gpu
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH -o /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/logs/04_chrombpnet_train_%j.out
#SBATCH -e /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/logs/04_chrombpnet_train_%j.err

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
BIAS_MODEL=$BASE/03_bias_model/models/bias.h5

OUTDIR=$BASE/04_chrombpnet_model
mkdir -p $OUTDIR

echo "[$(date)] === start training ChromBPNet ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# Check bias model
if [ ! -f "$BIAS_MODEL" ]; then
    echo "ERROR: bias model پیدا نشد: $BIAS_MODEL" >&2
    exit 1
fi
echo "Bias model: $BIAS_MODEL  ($(du -sh $BIAS_MODEL | cut -f1))"

# GPU test
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
    chrombpnet pipeline \
        -ibam $BAM \
        -d "ATAC" \
        -g $GENOME \
        -c $CHROMSIZES \
        -p $PEAKS \
        -n $NONPEAKS \
        -fl $FOLD \
        -b $BIAS_MODEL \
        -o $OUTDIR

echo "[$(date)] === ChromBPNet training is finished ==="

echo "--- output files ---"
for f in \
    "$OUTDIR/models/chrombpnet.h5" \
    "$OUTDIR/models/chrombpnet_nobias.h5" \
    "$OUTDIR/models/bias_model_scaled.h5" \
    "$OUTDIR/evaluation/overall_report.html"
do
    if [ -f "$f" ]; then
        echo "✓ $(basename $f)  ($(du -sh $f | cut -f1))"
    else
        echo "✗ is not found: $f"
    fi
done

echo "Genral report: $OUTDIR/evaluation/overall_report.html"
echo "loss performanc report:  cat $OUTDIR/logs/chrombpnet.log"
