# Go to project directory
cd /mnt/e/Projects/down_syndrome_PBMC

# Build HISAT2 index
hisat2-build results/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
results/reference/GRCh38_index

# Alignment SRR11856162
hisat2 -p 4 -x results/reference/GRCh38_index \
-1 data/fastq-trimmed/SRR11856162_1_paired.fastq \
-2 data/fastq-trimmed/SRR11856162_2_paired.fastq \
-S results/hisat2/SRR11856162.sam

# Alignment SRR11856163
hisat2 -p 4 -x results/reference/GRCh38_index \
-1 data/fastq-trimmed/SRR11856163_1_paired.fastq \
-2 data/fastq-trimmed/SRR11856163_2_paired.fastq \
-S results/hisat2/SRR11856163.sam

# SAM → BAM
samtools view -bS results/hisat2/SRR11856162.sam > results/hisat2/SRR11856162.bam
samtools view -bS results/hisat2/SRR11856163.sam > results/hisat2/SRR11856163.bam

# Sort BAM
samtools sort results/hisat2/SRR11856162.bam -o results/hisat2/SRR11856162.sorted.bam
samtools sort results/hisat2/SRR11856163.bam -o results/hisat2/SRR11856163.sorted.bam

# Index BAM
samtools index results/hisat2/SRR11856162.sorted.bam
samtools index results/hisat2/SRR11856163.sorted.bam

# FeatureCounts
featureCounts -T 4 -p -t exon -g gene_id \
-a results/reference/Homo_sapiens.GRCh38.109.gtf \
-o results/counts/gene_counts.txt \
results/hisat2/SRR11856162.sorted.bam \
results/hisat2/SRR11856163.sorted.bam
