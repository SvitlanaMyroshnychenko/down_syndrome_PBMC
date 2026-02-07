#!/bin/bash

# Check versions
fastqc --version
multiqc --version

# Create output directories
# mkdir -p results/fastqc-raw
# mkdir -p results/multiqc-raw

# Run FastQC on all raw FASTQ files
fastqc data/fastq-raw/*.fastq.gz -o results/fastqc-raw

# Generate MultiQC report
multiqc results/fastqc-raw -o results/multiqc-raw

# FastQC for trimmed reads
# fastqc data/fastq-trimmed/*.fastq.gz -o results/fastqc-trimmed

# FastQC for trimmed paired reads
fastqc data/fastq-trimmed/*_paired.fastq.gz -o results/fastqc-trimmed


# MultiQC for trimmed reads
multiqc results/fastqc-trimmed -o results/multiqc-trimmed

#### TRIMMOMATIC

java -jar E:\Programs\Trimmomatic\trimmomatic-0.40.jar PE -threads 4 ^
data\fastq-raw\SRR11856166_1.fastq.gz ^
data\fastq-raw\SRR11856166_2.fastq.gz ^
data\fastq-trimmed\SRR11856166_1_paired.fastq.gz ^
data\fastq-trimmed\SRR11856166_1_unpaired.fastq.gz ^
data\fastq-trimmed\SRR11856166_2_paired.fastq.gz ^
data\fastq-trimmed\SRR11856166_2_unpaired.fastq.gz ^
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 ^
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
