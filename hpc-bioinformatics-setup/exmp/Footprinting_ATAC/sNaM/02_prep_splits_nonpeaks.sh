cd /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted
mkdir -p 02_prep/splits

grep -E "^chr([1-9]|1[0-9]|X|Y)\b" \
    /mnt/archive/farhadie/ref/mm10/mm10.chrom.sizes \
    > 02_prep/mm10.main.chrom.sizes

cat 02_prep/mm10.main.chrom.sizes
