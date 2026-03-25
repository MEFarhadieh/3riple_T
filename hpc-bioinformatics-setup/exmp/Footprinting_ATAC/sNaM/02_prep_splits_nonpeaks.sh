cd /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted
mkdir -p 02_prep/splits

grep -E "^chr([1-9]|1[0-9]|X|Y)\b" \
    /mnt/archive/farhadie/ref/mm10/mm10.chrom.sizes \
    > 02_prep/mm10.main.chrom.sizes

cat 02_prep/mm10.main.chrom.sizes

# ─────────────────────────────────────────────
# build fold_0.json — chrombpnet 
#   test:  chr1, chr8   (~10%)
#   valid: chr2, chrX   (~10%)
#   train: rest         (~80%)
# ─────────────────────────────────────────────

singularity exec --bind /mnt/archive \
    /mnt/archive/farhadie/container/chrombpnet.sif \
    chrombpnet prep splits \
        -c /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/02_prep/mm10.main.chrom.sizes \
        -tcr chr1 chr8 \
        -vcr chr2 chrX \
        -op /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/02_prep/splits/fold_0


cat /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/02_prep/splits/fold_0.json



grep -E "^chr([1-9]|1[0-9]|X|Y)\b" \
    01_peaks/sNaM_peaks_no_blacklist.bed \
    > 01_peaks/sNaM_peaks_no_blacklist_mainchrom.bed

wc -l 01_peaks/sNaM_peaks_no_blacklist_mainchrom.bed

# ─────────────────────────────────────────────
# build nonpeaks — chrombpnet 
# ─────────────────────────────────────────────

singularity exec --bind /mnt/archive \
    /mnt/archive/farhadie/container/chrombpnet.sif \
    chrombpnet prep nonpeaks \
        -g /mnt/archive/farhadie/ref/mm10/mm10.fa \
        -p /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/01_peaks/sNaM_peaks_no_blacklist_mainchrom.bed \
        -c /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/02_prep/mm10.main.chrom.sizes \
        -fl /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/02_prep/splits/fold_0.json \
        -br /mnt/archive/farhadie/ref/mm10/mm10-blacklist.v2.bed \
        -npr 2 \
        -o /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted/02_prep/sNaM

wc -l 02_prep/sNaM_negatives.bed


grep -E "^chr([1-9]|1[0-9]|X|Y|M)\b" \
    /mnt/archive/farhadie/ref/mm10/mm10.chrom.sizes \
    > 02_prep/mm10.main.chrom.sizes

