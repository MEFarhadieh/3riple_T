# Prepare GWAS input for SMR analysis

cat > annotate_complete.sh << 'EOF'
#!/bin/bash

INPUT="II_GCTA_COJO.txt"
DBSNP="$LAB_ARCHIVE/ref/hg38/dbsnp_ref/GCF_000001405.25.gz"
OUTPUT="II_GCTA_COJO_with_rsids.txt"

echo "=== Step 1: Creating chromosome mapping ==="
cat > chrom_map.txt << 'MAPEND'
1 NC_000001.10
2 NC_000002.11
3 NC_000003.11
4 NC_000004.11
5 NC_000005.9
6 NC_000006.11
7 NC_000007.13
8 NC_000008.10
9 NC_000009.11
10 NC_000010.10
11 NC_000011.9
12 NC_000012.11
13 NC_000013.10
14 NC_000014.8
15 NC_000015.9
16 NC_000016.9
17 NC_000017.10
18 NC_000018.9
19 NC_000019.9
20 NC_000020.10
21 NC_000021.8
22 NC_000022.10
X NC_000023.10
Y NC_000024.9
MAPEND

echo "=== Step 2: Converting GWAS to VCF ==="
cat > temp.vcf << 'HEADER'
##fileformat=VCFv4.2
##INFO=<ID=BETA,Number=1,Type=Float,Description="Effect size">
##INFO=<ID=SE,Number=1,Type=Float,Description="Standard error">
##INFO=<ID=PVAL,Number=1,Type=Float,Description="P-value">
##INFO=<ID=FREQ,Number=1,Type=Float,Description="Allele frequency">
##INFO=<ID=N,Number=1,Type=Integer,Description="Sample size">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
HEADER

tail -n +2 ${INPUT} | awk -F'\t' '{
    split($1, a, ":");
    printf "%s\t%s\t.\t%s\t%s\t.\t.\tBETA=%s;SE=%s;PVAL=%s;FREQ=%s;N=%s\n", 
           a[1], a[2], a[3], a[4], $5, $6, $7, $4, $8
}' >> temp.vcf

echo "=== Step 3: Compressing and renaming chromosomes ==="
bgzip -f temp.vcf
tabix -f -p vcf temp.vcf.gz

bcftools annotate --rename-chrs chrom_map.txt temp.vcf.gz -Oz -o temp_renamed.vcf.gz
tabix -f -p vcf temp_renamed.vcf.gz

echo "=== Step 4: Annotating with rsIDs ==="
bcftools annotate \
    -a ${DBSNP} \
    -c ID \
    -Oz \
    -o annotated.vcf.gz \
    temp_renamed.vcf.gz

echo "=== Step 5: Converting back to table ==="
bcftools query \
    -f '%ID\t%REF\t%ALT\t%INFO/FREQ\t%INFO/BETA\t%INFO/SE\t%INFO/PVAL\t%INFO/N\n' \
    annotated.vcf.gz > temp_data.txt

echo -e "SNP\tA1\tA2\tfreq\tb\tse\tP\tN" > ${OUTPUT}
cat temp_data.txt >> ${OUTPUT}

echo "=== Step 6: Results ==="
total=$(tail -n +2 ${OUTPUT} | wc -l)
rsids=$(grep -c '^rs' ${OUTPUT} || echo 0)

echo "Total variants: $total"
echo "rsIDs assigned: $rsids"
echo "Output file: ${OUTPUT}"

# Cleanup
rm temp.vcf.gz temp.vcf.gz.tbi temp_renamed.vcf.gz temp_renamed.vcf.gz.tbi annotated.vcf.gz temp_data.txt chrom_map.txt

echo "DONE!"
EOF

chmod +x annotate_complete.sh
bash annotate_complete.sh
