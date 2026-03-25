#!/bin/bash
#SBATCH -J sNaM_peak_calling
#SBATCH -p hpc
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -o /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/logs/01_peak_calling_%j.out
#SBATCH -e /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/logs/01_peak_calling_%j.err

set -euo pipefail


export PYTHONPATH=/usr/people/EDVZ/farhadie/.local/lib/python3.10/site-packages:${PYTHONPATH:-}
export PATH=/usr/people/EDVZ/farhadie/.local/bin:$PATH

echo "Python: $(which python3)"
echo "macs2:  $(which macs2)"
python3 -c "import MACS2; print('MACS2 import OK')"

# ─────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────
BASE=/mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted
SIF=/mnt/archive/farhadie/container/chrombpnet.sif

BAM=$BASE/sNaM_sorted.bam
BLACKLIST=/mnt/archive/farhadie/ref/mm10/mm10-blacklist.v2.bed
CHROMSIZES=/mnt/archive/farhadie/ref/mm10/mm10.chrom.sizes

OUTDIR=$BASE/01_peaks
mkdir -p $OUTDIR $BASE/logs

echo "[$(date)] === شروع peak calling ==="

# ─────────────────────────────────────────────
# peak calling with macs2
# ─────────────────────────────────────────────
macs2 callpeak \
    -t $BAM \
    -f BAM \
    -n sNaM \
    --outdir $OUTDIR/macs2_out \
    --nomodel \
    --shift -100 \
    --extsize 200 \
    -p 0.01 \
    --keep-dup all \
    -g mm \
    --call-summits

echo "[$(date)] MACS2 is finished"
echo "total raw peaks : $(wc -l < $OUTDIR/macs2_out/sNaM_peaks.narrowPeak)"

# ─────────────────────────────────────────────
#  blacklist filitration with bedtools 
# ─────────────────────────────────────────────
bedtools slop \
    -i $BLACKLIST \
    -g $CHROMSIZES \
    -b 1057 \
    > $OUTDIR/blacklist_slop1057.bed

bedtools intersect \
    -v \
    -a $OUTDIR/macs2_out/sNaM_peaks.narrowPeak \
    -b $OUTDIR/blacklist_slop1057.bed \
    > $OUTDIR/sNaM_peaks_no_blacklist.bed

echo "[$(date)] filtration of blacklist is finished"
echo "total peaks post blacklist filitration: $(wc -l < $OUTDIR/sNaM_peaks_no_blacklist.bed)"

echo "[$(date)] === Step 1 is finished ==="
echo "output: $OUTDIR/sNaM_peaks_no_blacklist.bed"
