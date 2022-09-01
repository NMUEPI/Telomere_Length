# Telomere Length Analysis Based on Whole-genome Sequecing

This repository contains custom scripts used for whole-genome sequencing processing, telomere length estimation, single-cell multiomics data analysis, together with transcription factor analysis.

## Overview

### Whole-genome Sequencing 
Mapping read sequences to the human reference genome.

### Telomere Length Estimation
Telseq *(https://github.com/zd1/telseq)* and Qmotif *(https://github.com/AdamaJava/adamajava/tree/master/qmotif)* were applied to estimate the telomere length based on the whole-genome sequencing data.

### Single-cell Multiomics Analysis
Datasets can be accessed in the NCBI Gene Expression Omnibus (GEO) repository under their GSE numbers: accession numbers of mouse embryo scRNA-seq data are GSE136714; accession numbers of mouse embryo scNOME-seq data are GSE136715.

### Singel-cell Transcription Factor Analysis
Based on the single-cell RNA-seq data, R package named Scenic *(https://github.com/aertslab/SCENIC)* was used to perform transcript factor analysis.
