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
