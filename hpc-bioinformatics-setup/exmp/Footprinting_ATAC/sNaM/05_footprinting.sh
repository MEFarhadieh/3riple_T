#!/bin/bash
#SBATCH -J sNaM_footprint
#SBATCH -p hpc
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -o /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/logs/05_footprinting_%j.out
#SBATCH -e /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/logs/05_footprinting_%j.err

set -euo pipefail

# ─────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────
BASE=/mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted
SIF=/mnt/archive/farhadie/container/chrombpnet.sif

MODEL=$BASE/04_chrombpnet_model/models/chrombpnet_nobias.h5
NONPEAKS=$BASE/02_prep/sNaM_negatives.bed
GENOME=/mnt/archive/farhadie/ref/mm10/mm10.fa
FOLD=$BASE/02_prep/splits/fold_0.json
MOTIFS_TSV=$BASE/05_footprinting/motifs_to_pwm.tsv

OUTDIR=$BASE/05_footprinting
mkdir -p $OUTDIR

echo "[$(date)] === start marginal footprinting ==="

# ─────────────────────────────────────────────
# Build motifs TSV
# ─────────────────────────────────────────────
cat > $MOTIFS_TSV << 'MOTIF_EOF'
RUNX3	TGTGGT
SPI1	GGAAGTG
IRF4	GAGGAAGTG
CEBPB	ATTGCGCAAT
AP1	TGASTCAGC
SMAD3	GTCTGNNCAGAC
SMAD4	GTCTGNNCAGAC
SMAD2	AGTATGTCTAGC
MOTIF_EOF

echo "motifs file:"
cat $MOTIFS_TSV

# ─────────────────────────────────────────────
# Marginal Footprinting
# ─────────────────────────────────────────────
singularity exec --bind /mnt/archive --bind /cl_tmp $SIF \
    chrombpnet footprints \
        -m $MODEL \
        -r $NONPEAKS \
        -g $GENOME \
        -fl $FOLD \
        -op $OUTDIR/sNaM \
        -pwm_f $MOTIFS_TSV

echo "[$(date)] === Footprinting is finished==="
echo "Output directory:"
ls -lh $OUTDIR/
