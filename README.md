# MendelPeaG2P

## Project Overview
This repository contains all code and analysis pipelines related to our research published in **Biorxiv**.

## Directory Structure and Functions

### 00.SNP_calling_and_QC
This directory contains a complete pipeline from sequencing data to high-quality SNP sets, including:
- Raw sequencing data quality control
- Alignment to reference genome
- SNP detection and genome variation annotation
- Multi-step SNP quality filtering (missing rate, frequency, linkage disequilibrium pruning, etc.)

### 01.HAPMAP_pipeline
This directory contains workflows for constructing haplotype maps, used to analyze the genetic diversity and evolutionary relationships of pea germplasm resources.

### 02.PCA_t-SNE_and_population_structure
This directory contains scripts for population structure analysis, including:
- Principal Component Analysis (PCA)
- t-SNE dimensionality reduction analysis
- Population structure analysis methods such as ADMIXTURE

### 03.GWAS_pipeline
This directory contains a complete pipeline for Genome-Wide Association Studies (GWAS), including:
- Main GWAS analysis scripts
- Manhattan plot and QQ plot drawing functions
- Implementation of various statistical models

### 04.HAPPE_for_Pea
This directory contains the implementation of HAPPE software specifically for pea, including:
- Haplotype construction and visualization
- Functional variation annotation and filtering
- Hierarchical clustering analysis
- Data visualization and result export

### 05.CNV_pipeline
This directory contains pipelines for Copy Number Variation (CNV) detection and analysis.

### 06.haplotype-GWAS
This directory contains GWAS analysis methods based on haplotypes, using haplotypes as markers for association analysis.

### 07.BSA.analysis
This directory contains Bulked Segregant Analysis (BSA) pipelines for pea, used to rapidly locate genomic regions associated with target traits.

## Usage
The scripts in each directory should be executed in numerical order. Please refer to the script comments in each directory for detailed parameter descriptions and usage instructions. Most scripts need to be run in a Linux/Unix environment and depend on common bioinformatics tools.

## Citation
Feng, C., Chen, B., Hofer, J., Shi, Y., Jiang, M., Song, B., ... & Cheng, S. (2024). Genomic and Genetic Insights into Mendel’s Pea Genes. bioRxiv, 2024-05.

## Contact
For any questions or suggestions, please contact us at:
[fengcong@caas.cn]

## License
This project is licensed under the MIT License. See the LICENSE file for details.
