# Comparative Meta-Analysis of Abiotic Stress Responses in Maize (*Zea mays* L.)

This repository contains the full computational pipeline and analytical scripts used for my Bachelor's Thesis titled *"Differential Gene Expression Analysis of
Maize Abiotic Stress Response Using Public RNA-Seq Datasets"*. 

The project performs a transcriptomic meta-analysis to identify conserved gene expression patterns in maize under drought, cold, and heat stress, and further evaluates their phenotypic relevance using GWAS data.

## Repository Structure

The repository is organized to follow the logical flow of the thesis:

- `/Scripts/`: All R scripts for data processing and statistical analysis.
- `/Metadata/`: Sample annotations, GEO accessions
- `/Results/`: Processed DEG lists and enrichment results.

## Analysis Workflow

### 1. Data Retrieval & Preprocessing
Raw RNA-seq data (FASTQ) were retrieved from NCBI GEO/SRA using the **SRA Toolkit**.
- **Pipeline:** Preprocessing and quantification were performed using the https://github.com/nf-core/rnaseq pipeline (Nextflow).

### 2. Differential Expression Analysis
Statistical testing was performed using **DESeq2**.
- **Script:** `Scripts/01_DESeq2_Analysis.R`
- **Key Steps:** Median-of-ratios normalization, `apeglm` LFC shrinkage, and Benjamini-Hochberg FDR correction.

### 3. Functional Enrichment
Biological interpretation of DEGs using GO and KEGG.
- **Script:** `Scripts/02_GO_Enrichment.R`and `Scripts/03_KEGG_Enrichment.R`
- **Tools:** `gprofiler2` for over-representation testing.

### 4. Cross-Study & Cross-Stress Meta-Analysis
Identifying the 496-gene conserved "core" between drought and cold.
- **Script:** `Scripts/04_Meta_Analysis_Drought.R`, `Scripts/05_Meta_Analysis_Cold.R` ,`Scripts/06_Meta_Analysis_Heat.R`,`Scripts/07_Enrichment_Overlapping_Genes.R` ,`Scripts/08_Functional_Meta_Analysis.R`
- **Visuals:** Venn diagrams, UpSet plots, and functional similarity heatmaps.

### 5. GWAS Integration
Mapping stress genes to phenotypic yield data using NAM population statistics.
- **Script:** `Scripts/09_Mapping_Phenotype.R`
- **Method:** Genomic interval intersection ($\pm 50$ kb windows) and Fisher's Exact Test.
