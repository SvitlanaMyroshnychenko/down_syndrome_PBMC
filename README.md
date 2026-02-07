## Down Syndrome PBMC RNA-seq Analysis

### Goal

To explore transcriptomic differences between Down Syndrome and healthy PBMC samples using RNA-seq data and a reproducible bioinformatics workflow.

### Note

> This project uses publicly available datasets for learning and demonstration purposes only.

------------------------------------------------------------------------

### Data

**Source:** [NCBI GEO – GSE151282](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151282)\
Raw sequencing data available via SRA (SRP264957).

**Project:** GSE151282 / SRP264957

**Data type:** RNA-seq, paired-end, PBMC samples

**Selected samples (4 total):**

| SRR         | Condition     | Sex    |
|-------------|---------------|--------|
| SRR11856162 | Down syndrome | Male   |
| SRR11856166 | Down syndrome | Female |
| SRR11856164 | Healthy       | Male   |
| SRR11856165 | Healthy       | Female |

**Format:** `.fastq.gz` (paired-end, compressed)

------------------------------------------------------------------------

### Raw files

`data/fastq-raw/`

SRR11856162_1.fastq.gz\
SRR11856162_2.fastq.gz

SRR11856164_1.fastq.gz\
SRR11856164_2.fastq.gz

SRR11856165_1.fastq.gz\
SRR11856165_2.fastq.gz

SRR11856166_1.fastq.gz\
SRR11856166_2.fastq.gz

### Current Status

-   Downloaded raw FASTQ files for four samples:

    -   `SRR11856162` (Down syndrome)
    -   `SRR11856166` (Down syndrome)
    -   `SRR11856164` (Healthy)
    -   `SRR11856165` (Healthy)

-   Performed quality control on raw reads using FastQC

-   Generated MultiQC summary report for raw reads\
    → `results/multiqc-raw/multiqc_report.html`

-   Performed trimming using Trimmomatic (quality + adapter trimming)

-   Generated FastQC and MultiQC reports for trimmed reads\
    → `results/multiqc-trimmed/multiqc_report.html`

-   After comparison of raw vs trimmed reports, trimming did not significaly improve quality metrics.\
    That's why raw reads were selected for downstream alignment and analysis.

------------------------------------------------------------------------

### Quality Control

FastQC and MultiQC reports were generated for both raw and trimmed reads.

**Reports locations:**

-   FastQC: `results/fastqc-raw/` (raw), `results/fastqc-trimmed/` (trimmed)\
-   MultiQC: `results/multiqc-raw/multiqc_report.html` (raw), `results/multiqc-trimmed/multiqc_report.html` (trimmed)

**Observations:**

-   **Raw reads:**
    -   Overall per-base quality is high across all samples (green in Per Sequence Quality Scores)
    -   No significant adapter contamination detected (green in Adapter Content)
    -   GC content deviates slightly from expected (red in Per Sequence GC Content)
    -   High sequence duplication levels observed (red in Sequence Duplication Levels)
    -   Sequence lengths are uniform at 125 bp (green in Sequence Length Distribution)
-   **Trimmed reads:**
    -   Low-quality bases from read ends were removed using Trimmomatic
    -   Total read counts reduced by \~5–10% due to trimming
    -   Sequence lengths became variable after trimming (expected outcome)
    -   Overall, trimming did not significaly improve quality metrics

**Conclusion:** FastQC and MultiQC reports are available for both raw and trimmed reads. For downstream alignment, quantification and all future analyses, raw paired-end reads will be used.

### Alignment (HISAT2)

-   Downloaded GRCh38 reference genome (Ensembl release 109)

-   Downloaded corresponding GTF annotation file

-   Built HISAT2 index\
    → `results/reference/GRCh38_index.*.ht2`

-   Aligned raw paired-end reads to reference genome

-   Converted SAM → BAM → sorted BAM

-   Indexed BAM files

Final alignment files:\
→ `results/hisat2/*_sorted.bam`\
→ `results/hisat2/*_sorted.bam.bai`

------------------------------------------------------------------------

### Quantification (featureCounts)

-   Generated gene-level count matrix using exon-based counting
-   Counted paired-end reads (`--countReadPairs`)
-   Annotation: `Homo_sapiens.GRCh38.109.gtf`

Output files: - Gene count matrix → `results/featurecounts_counts.txt` - Assignment summary → `results/featurecounts_counts.txt.summary`

------------------------------------------------------------------------

### Next Steps

-   Variance stabilizing transformation (VST)
-   PCA for sample clustering
-   Differential expression analysis using DESeq2
