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

`results/fastqc/`

SRR11856162_1_fastqc.html

SRR11856162_2_fastqc.html

SRR11856163_1_fastqc.html

SRR11856163_2_fastqc.html

### Current Status

-   Downloaded raw fastq files for SRR11856162 and SRR11856163
-   Performed quality control using FastQC
-   Results (HTML reports and ZIP files) are saved in `results/fastqc/`
-   MultiQC summary report compiled for all samples.

### Quality Control

Observations:

-   Overall per-base quality is high across all samples (green in Per Sequence Quality Scores).

-   No significant adapter contamination detected (green in Adapter Content).

-   GC content deviates slightly from expected (red in Per Sequence GC Content).

-   High sequence duplication levels observed (red in Sequence Duplication Levels).

-   Sequence lengths are as expected (green in Sequence Length Distribution).

-   Reports are available in `results/multiqc/multiqc_report.html`.

### Next Steps

-   Trimming low-quality reads
-   Alignment to reference genome
-   Count matrix generation
-   Exploratory analysis (VST, PCA)
-   Differential expression analysis (DESeq2)
-   Clustering and visualization of top DE genes
