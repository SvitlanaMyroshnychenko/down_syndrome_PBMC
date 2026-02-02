## Down Syndrome PBMC RNA-seq Analysis

### Goal

The goal of this project is to explore transcriptomic features of PBMC cells from Down Syndrome patients, focusing on differentially expressed genes and cell clustering to understand biological pathways involved in the disease.

### Data

**Source:** NCBI GEO / SRA

**Project:** GSE151282 / SRP264957

**Data type:** RNA-seq, paired-end, PBMC samples

**Downloaded:** first 2 samples for pilot analysis

**Format:** .fastq raw sequencing files

**FastQC:** quality control reports in .html

#### Raw files

`data/fastq-raw/`

SRR11856162_1.fastq

SRR11856162_2.fastq

SRR11856163_1.fastq

SRR11856163_2.fastq

### Analysis Workflow

-   Raw data download from SRA

-   Quality control using FastQC

-   Adapter and quality trimming

-   Post-trimming quality control

-   Summary report using MultiQC

-   (Planned) Alignment and quantification

-   (Planned) Differential expression analysis using DESeq2

### Current Status

-   Downloaded raw FASTQ files for `SRR11856162` and `SRR11856163`
-   Performed initial quality control on raw reads using FastQC
-   Trimmed low-quality bases and adapter sequences from all reads
-   Performed FastQC on trimmed reads to assess post-trimming quality
-   Compiled MultiQC summary reports for both raw and trimmed reads
-   All results (HTML reports and ZIP files) are saved in `results/fastqc/` and `results/multiqc-trimmed/`

### Quality Control

#### Reports

FastQC reports (raw and trimmed):

`results/fastqc/`

`results/fastqc-trimmed`

MultiQC summary reports:

`results/multiqc/multiqc_report.html`

`results/multiqc-trimmed/multiqc_report.html`

### Observations

#### Raw reads

-   Overall per-base quality is high across all samples (green in Per Sequence Quality Scores).

-   No significant adapter contamination detected (green in Adapter Content).

-   GC content deviates slightly from expected (red in Per Sequence GC Content).

-   High sequence duplication levels observed (red in Sequence Duplication Levels).

-   Sequence lengths are as expected (green in Sequence Length Distribution).

#### Trimmed reads

-   Removal of low-quality bases from 3′ ends

-   Variable read length after trimming (expected outcome)

-   Reduction in total read count (\~9%)

-   Improved read quality suitable for downstream alignment and quantification

-   Based on these observations, trimmed reads were selected for downstream analysis.

### Next Steps

-   Alignment to reference genome
-   Count matrix generation
-   Exploratory analysis (VST, PCA)
-   Differential expression analysis (DESeq2)
-   Clustering and visualization of top DE genes

#### Notes

This project is intended as a learning and demonstration pipeline and is not designed to reproduce the full-scale analysis of the original study.
