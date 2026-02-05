#!/bin/bash
# Download SRA files

cd /mnt/e/Projects/down_syndrome_PBMC

# download SRA files from the assession list
while read srr; do
    fastq-dump --split-files --gzip $srr -O data/fastq-raw
done < data/fastq-raw/SRR_Acc_List.txt

# download a specific SRA file from the assession list (optional)
fastq-dump --split-files --gzip SRR11856166 -O data/fastq-raw
