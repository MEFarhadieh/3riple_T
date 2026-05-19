singularity exec --bind /mnt/archive \
    /mnt/archive/farhadie/container/chrombpnet.sif \
    modisco report \
        -i 18_modisco_1M/modisco_results_1M.h5 \
        -o 18_modisco_1M/report/ \
        -s 18_modisco_1M/report/ \
        -m ../sNaM_sorted/18_modisco_1M/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt
