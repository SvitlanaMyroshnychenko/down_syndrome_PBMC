#!/bin/bash

# Go to project directory
cd /mnt/e/Projects/down_syndrome_PBMC


# 1. REFERENCE PREPARATION

# Create reference folder
mkdir -p results/reference
cd results/reference

# Download genome (FASTA)
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Download GTF annotation
wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

# Unzip files
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.109.gtf.gz


# 2. BUILD HISAT2 INDEX

hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38_index


# 3. ALIGNMENT

cd /mnt/e/Projects/down_syndrome_PBMC
mkdir -p results/hisat2

for fwd in data/fastq-raw/*_1.fastq.gz; do
    base=$(basename "$fwd" _1.fastq.gz)
    rev="data/fastq-raw/${base}_2.fastq.gz"

    hisat2 -p 4 -x results/reference/GRCh38_index \
        -1 "$fwd" -2 "$rev" | \
    samtools view -@ 4 -bS - | \
    samtools sort -@ 4 -o "results/hisat2/${base}_sorted.bam"

    samtools index "results/hisat2/${base}_sorted.bam"
done


# 4. FEATURECOUNTS

featureCounts -T 4 -p --countReadPairs -t exon -g gene_id \
-a results/reference/Homo_sapiens.GRCh38.109.gtf \
-o results/featurecounts_counts.txt \
results/hisat2/*_sorted.bam
