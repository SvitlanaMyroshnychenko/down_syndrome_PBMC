# Down Syndrome PBMC RNA-seq Analysis

### Project Goal

The primary objective of this project was to implement an end-to-end RNA-seq bioinformatics pipeline to identify transcriptomic signatures in Peripheral Blood Mononuclear Cells (PBMC) of individuals with **Down Syndrome**.

------------------------------------------------------------------------

### Dataset Summary

-   **Source:** [NCBI GEO – GSE151282](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151282)
-   **Design:** 4 paired-end samples, balanced by condition and sex.
-   **Organism:** *Homo sapiens* (GRCh38.109)

| Sample ID   | Condition     | Sex    |
|-------------|---------------|--------|
| SRR11856162 | Down Syndrome | Male   |
| SRR11856166 | Down Syndrome | Female |
| SRR11856164 | Healthy       | Male   |
| SRR11856165 | Healthy       | Female |

------------------------------------------------------------------------

### Bioinformatics Pipeline

1.  **Quality Control:** `FastQC` & `MultiQC`. Raw reads showed high quality, eliminating the need for aggressive trimming.
2.  **Alignment:** `HISAT2`. Splice-aware mapping to the GRCh38 human primary assembly.
3.  **Quantification:** `featureCounts` (Subread package). Gene-level counting based on Ensembl GTF annotation.
4.  **Differential Expression:** `DESeq2` (R-based). Statistical modeling, VST transformation and p-value adjustment.

------------------------------------------------------------------------

### Key Findings

-   **Transcriptomic Separation.** PCA analysis (PC1: 55.5% variance) and Hierarchical Clustering show a clear, distinct separation between Down Syndrome and Healthy groups.
-   **Chromosome 21 Validation.** Significant up-regulation of the MX1 gene was detected. Being located on Chromosome 21, its over-expression serves as a direct internal validation of the Trisomy 21 gene dosage effect.
-   **Biomarkers.** Identified high-confidence Differentially Expressed Genes (DEGs) including NRG1 (neurodevelopmental signaling) and HBG1 (hemoglobin subunit gamma-1).

------------------------------------------------------------------------
