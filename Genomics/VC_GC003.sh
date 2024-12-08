#!/bin/bash

bwa mem hg19.fasta sample_1.fastq sample_2.fastq > sample.sam

samtools view -Sb sample.sam | samtools sort -o sample_sorted.bam

gatk MarkDuplicates \
   -I sample_sorted.bam \
   -O sample_dedup.bam \
   -M sample_metrics.txt


gatk BaseRecalibrator \
   -I sample_dedup.bam \
   -R hg19.fasta \
   --known-sites dbsnp.vcf \
   -O recal_data.table
   
gatk ApplyBQSR \
    -R hg19.fasta \
    -I sample_dedup.bam \
    --bqsr-recal-file recal_data.table \
    -O sample_recal.bam

gatk VariantFiltration \
    -R hg19.fasta \
    -V sample_raw_variants.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0" \
    --filter-name "low_quality" \
    -O sample_filtered_variants.vcf
