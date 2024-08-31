################################################################################

# ChIP-seq analysis of PRDM16 KO/WT in cortical neurons

################################################################################

# 1- Planning and organization

################################################################################

cd ~
mkdir chipseq
cd chipseq
mkdir raw_data reference_data scripts logs meta results
tree

# Download raw data from ENA
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP134733&result=read_run&fields=run_accession,fastq_ftp&format=tsv&download=true" > ena_fastq_urls.tsv
cut -f 2 ena_fastq_urls.tsv | tr ";" "\n" | while read url; do wget -P ~/chipseq/raw_data "$url"; done
gunzip ~/chipseq/raw_data/*.fastq.gz

# Use https://www.ncbi.nlm.nih.gov/sra?term=SRP134733 to geratr
while IFS=$'\t' read -r old_name new_name; do
    mv ~/chipseq/raw_data/"$old_name" ~/chipseq/raw_data/"$new_name"
done < rename_list.tsv

################################################################################

# 2- Quality control of sequence reads

################################################################################

mkdir -p ~/chipseq/results/fastqc
fastqc -o ~/chipseq/results/fastqc/ -t 8 ~/chipseq/raw_data/*chip.fastq.gz

################################################################################

# 3- Alignment to Genome

################################################################################

# Build index
cd ~/chipseq/reference_data/
bowtie2-build ~/chipseq/reference_data/mm10_reference_genome.fa mm10


input_dir=~/chipseq/raw_data/
output_dir_sam=~/chipseq/results/bowtie2/
output_dir_log=~/chipseq/logs/
index_path=~/chipseq/reference_data/mm10

for fastq_file in "$input_dir"*chip.fastq.gz; do
  base_name=$(basename "$fastq_file" .fastq.gz)
  sam_output="${output_dir_sam}${base_name}.sam"
  log_output="${output_dir_log}${base_name}_bowtie2.log"
  # Run aligner
  bowtie2 -p 8 -q --local \
        -x "$index_path" \
        -U "$fastq_file" \
        -S "$sam_output" 2> "$log_output"
done

# Convert SAM to BAM
input_dir=~/chipseq/results/bowtie2/
output_dir_bam=~/chipseq/results/bowtie2/

for sam_file in "$input_dir"*.sam; do
  base_name=$(basename "$sam_file" .sam)
  bam_output="${output_dir_bam}${base_name}.bam"
  # Run samtools
  samtools view -h -S -b \
  -o "$bam_output" \
  "$sam_file"
done

################################################################################

# 4- Filtering BAM files

################################################################################

# Sort BAM files by genomic coordinates
input_dir=~/chipseq/results/bowtie2/
output_dir_bam=~/chipseq/results/bowtie2/

for bam_file in "$input_dir"*.bam; do
  base_name=$(basename "$bam_file" .bam)
  bam_output="${output_dir_bam}${base_name}_sorted.bam"
  samtools sort "$bam_file" -o "$bam_output"
done

# Filter the reads to keep only uniquely mapping reads
for bam_file in "$input_dir"*sorted.bam; do
  base_name=$(basename "$bam_file" _sorted.bam)
  bam_output="${output_dir_bam}${base_name}_dedup.bam"
  sambamba view -h -t 8 -f bam \
  -F "[XS] == null and not unmapped and not duplicate" \
  "$bam_file" > "$bam_output"
done

# Filtering out Blacklisted Regions
for bam_file in "$input_dir"*dedup.bam; do
  base_name=$(basename "$bam_file" _dedup.bam)
  bam_output="${output_dir_bam}${base_name}_final.bam"
  bedtools intersect -v -abam "$bam_file" -b \
  mm10-blacklist.v2.bed > "$bam_output"
done

################################################################################

# 5- Peak calling

################################################################################

mkdir -p ~/chipseq/results/macs2
cd ~/chipseq/results/

macs2 callpeak -t /bowtie2/wt_sample1_chip_final.bam \
    -c /wt_sample1_input_final.bam \
    -f BAM -g mm \
    -n wt_sample1 \
    --outdir macs2 2> macs2/wt_sample1_macs2.log

macs2 callpeak -t /bowtie2/wt_sample2_chip_final.bam \
    -c /bowtie2/wt_sample2_input_final.bam \
    -f BAM -g mm \
    -n wt_sample2 \
    --outdir macs2 2> macs2/wt_sample2_macs2.log

mv macs2/*.log ../logs/

################################################################################

# 6- Handling peak calls

################################################################################

cd ~/chipseq/results/
# Filtering peaks overlapping with blacklist regions
bedtools intersect \
-v \
-a macs2/wt_sample1_peaks.narrowPeak \
-b ../reference_data/mm10-blacklist.v2.bed \
> macs2/wt_sample1_peaks_filtered.bed

bedtools intersect \
-v \
-a macs2/wt_sample2_peaks.narrowPeak \
-b ../reference_data/mm10-blacklist.v2.bed \
> macs2/wt_sample2_peaks_filtered.bed

# Finding overlapping peaks between replicates
bedtools intersect \
-wo -f 0.3 -r \
-a macs2/wt_sample1_peaks_filtered.bed \
-b macs2/wt_sample2_peaks_filtered.bed \
> macs2/wt_peaks_final.bed

wc -l ~/chipseq_workshop/results/macs2/wt_peaks_final.bed

################################################################################

# 7- Peak visualization

################################################################################

cd ~/chipseq/results/
mkdir -p visualization/bigWig

samtools index ~/chipseq/results/bowtie2/wt_sample2_chip_final.bam

bamCoverage -b ~/chipseq_workshop/results/bowtie2/wt_sample2_chip_final.bam \
-o ~/chipseq_workshop/results/visualization/bigWig/wt_sample2_chip.bw \
--binSize 20

computeMatrix reference-point --referencePoint center \
-b 4000 -a 4000 \
-R ~/chipseq_workshop/results/macs2/wt_peaks_final.bed \
-S visualization/bigWig/wt_sample1_chip.bw visualization/bigWig/wt_sample2_chip.bw \
--skipZeros \
-o ~/chipseq_workshop/results/visualization/wt_matrix.gz \
-p 8

mkdir ~/chipseq_workshop/results/visualization/figures
plotProfile -m ~/chipseq_workshop/results/visualization/wt_matrix.gz \
-out ~/chipseq_workshop/results/visualization/figures/plot1_wt_replicates.png \
--regionsLabel "" \
--perGroup \
--colors red blue \
--samplesLabel "WT_replicate1" "WT_replicate2" \
--refPointLabel "PRDM16 binding sites"

computeMatrix reference-point --referencePoint center \
  -b 4000 -a 4000 \
  -R ~/chipseq_workshop/results/macs2/wt_peaks_final.bed \
  -S visualization/bigWig/wt_sample2_chip.bw visualization/bigWig/ko_sample2_chip.bw \
  --skipZeros \
  -o visualization/wt_ko_matrix.gz \
  -p 6

plotProfile -m ~/chipseq_workshop/results/visualization/wt_ko_matrix.gz \
  -out ~/chipseq_workshop/results/visualization/figures/plot2_wt_ko.png \
  --regionsLabel "" \
  --perGroup \
  --colors blue red \
  --samplesLabel "WT" "KO" \
  --refPointLabel "PRDM16 binding sites"

plotProfile -m ~/chipseq_workshop/results/visualization/wt_matrix_allGenes_TSS.gz \
  -out ~/chipseq_workshop/results/visualization/figures/plot1_wt_TSS.png \
  --regionsLabel "" \
  --perGroup \
  --colors red blue \
  --samplesLabel "WT_replicate1" "WT_replicate2" \
  --refPointLabel "TSS" \
  --yMax 12

# Histone modifications and enhancers
plotProfile -m ~/chipseq_workshop/results/visualization/wt_encode_matrix.gz \
  -out ~/chipseq_workshop/results/visualization/figures/plot2_wt_encode.png \
  --regionsLabel "" \
  --perGroup \
  --colors blue green red orange \
  --samplesLabel "WT_replicate2" "H3K4me" "H3K27me3" "H3K27ac" \
  --refPointLabel "PRDM16 binding sites"
