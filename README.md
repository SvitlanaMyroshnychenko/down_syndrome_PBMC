## Down Syndrome PBMC RNA-seq Analysis

### Goal

#### Note

> This project uses public datasets for learning and demonstration purposes only.

### Data

**Source:** [NCBI GEO / SRA](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151282)\
Raw sequencing data available via SRA (SRP264957).

**Project:** GSE151282 / SRP264957

**Data type:** RNA-seq, paired-end, PBMC samples

**Downloaded:** first 2 samples for pilot analysis

**Format:** .fastq raw sequencing files

#### Raw files

`data/fastq-raw/`

SRR11856162_1.fastq

SRR11856162_2.fastq

SRR11856163_1.fastq

SRR11856163_2.fastq

### Current Status

-   Downloaded raw FASTQ files for `SRR11856162` and `SRR11856163`

-   Performed initial quality control on raw reads using FastQC

-   Trimmed low-quality bases and adapter sequences from all reads

-   Performed FastQC on trimmed reads to assess post-trimming quality

-   Compiled MultiQC summary reports for both raw and trimmed reads

-   All results (HTML reports and ZIP files) are saved in `results/fastqc/` and `results/multiqc-trimmed/`

-   Alignment to reference genome (HISAT2)

    -- Built HISAT2 index for GRCh38 → `results/reference/GRCh38_index.\*.ht2`

    -- Aligned trimmed reads → SAM files in `results/hisat2/`

    -- Converted SAM → sorted BAM → indexed BAM

-   Quantification (featureCounts)

    -- Generated gene-level count matrix → `results/counts/gene_counts.txt`

    -- Summary of read assignment → `results/counts/gene_counts_summary.txt`

-   Differential expression analysis using DESeq2 (planned)

    -- Count matrix ready for DESeq2 input

-   Next: exploratory analysis (PCA, VST) and identification of DE genes

### Quality Control

#### Reports

FastQC reports (raw and trimmed):

`results/fastqc/`

`results/fastqc-trimmed`

MultiQC summary reports:

`results/multiqc/multiqc_report.html`

`results/multiqc-trimmed/multiqc_report.html`

#### Quality Control Observations

##### Raw reads

-   Overall per-base quality is high across all samples (green in Per Sequence Quality Scores).

-   No significant adapter contamination detected (green in Adapter Content).

-   GC content deviates slightly from expected (red in Per Sequence GC Content).

-   High sequence duplication levels observed (red in Sequence Duplication Levels).

-   Sequence lengths are as expected (green in Sequence Length Distribution).

##### Trimmed reads

-   Removal of low-quality bases from 3′ ends

-   Variable read length after trimming (expected outcome)

-   Reduction in total read count (\~9%)

-   Improved read quality suitable for downstream alignment and quantification

-   Based on these observations, trimmed reads were selected for downstream analysis.


### Count Matrix Notes

Count matrix ready for differential expression analysis with DESeq2. Example summary (from featureCounts):

| Status | results/hisat2/SRR11856162.sorted.bam | results/hisat2/SRR11856163.sorted.bam |
|------------------------|------------------------|------------------------|
| Assigned | 26,750,659 | 17,076,077 |
| Unassigned_Unmapped | 1,523,505 | 1,470,596 |
| Unassigned_Read_Type | 0 | 0 |
| Unassigned_Singleton | 0 | 0 |
| Unassigned_MappingQuality | 0 | 0 |
| Unassigned_Chimera | 0 | 0 |
| Unassigned_FragmentLength | 0 | 0 |
| Unassigned_Duplicate | 0 | 0 |
| Unassigned_MultiMapping | 21,685,046 | 13,248,324 |
| Unassigned_Secondary | 0 | 0 |
| Unassigned_NonSplit | 0 | 0 |
| Unassigned_NoFeatures | 5,847,215 | 3,360,593 |
| Unassigned_Overlapping_Length | 0 | 0 |
| Unassigned_Ambiguity | 7,227,260 | 5,343,303 |

### Next Steps

-   Perform DESeq2 analysis to identify differentially expressed genes

-   Exploratory analysis: PCA, clustering, VST transformation
