# RNA-seq
Analysis of RNA-seq data of cassava (Manihot esculenta) plants infected with Xanthomonas phaseoli pv. manihotis from Zárate-Chaves et al. (2021).
This analysis is based on the tutorial by the Cebola Lab (https://github.com/CebolaLab/RNA-seq).

## Introduction

The first part of this analysis must be run in the terminal, by first creating an appropriate environment. In this case, we will create an environment called RNA-seq and install the required programs using Anaconda:

```
#Create environment
conda create -N RNA-seq

#Install required programs
conda install -n RNA-seq -c bioconda fastqc
conda install -n RNA-seq -c bioconda fastp
conda install -n RNA-seq -c bioconda multiqc
conda install -n RNA-seq -c bioconda star
conda install -n RNA-seq -c bioconda samtools
conda install -n RNA-seq -c bioconda deeptools
conda install -n RNA-seq -c bioconda salmon
conda install -n RNA-seq -c bioconda qualimap
conda install -n RNA-seq -c bioconda gffread

#For differential expression using DESeq2
conda create --name DEseq2 r-essentials r-base

conda install -n DEseq2 -c bioconda bioconductor-deseq2
conda install -n DEseq2 -c bioconda bioconductor-tximport 
conda install -n DEseq2 -c r r-ggplot2 

## Pre-alignment QC

## Post-alignment QC

![Genome Fraction Coverage](https://github.com/user-attachments/assets/4f5d8457-3527-4b1a-bf4f-bc613ae5fcdd)

![Insert Size Histogram](https://github.com/user-attachments/assets/a15d8fe3-b7a7-40fa-9cec-83adf0d13572)

## Visualization 

GC-bias.

Control

![SRR13313987 biasPlot](https://github.com/user-attachments/assets/851f67de-4733-4109-957e-24003bce5bc8)

UA681

![SRR13313993 biasPlot](https://github.com/user-attachments/assets/04aad013-3540-4b65-8083-f45b8e329867)

UA1061

![SRR13313990 biasPlot](https://github.com/user-attachments/assets/47cfec77-48c1-45e8-8e84-05ee68f221b2)

## Quantification

## Differential expression

![PCA](https://github.com/user-attachments/assets/6fc0c530-924b-452c-a8d7-0a040d751789)

![MA_UA681](https://github.com/user-attachments/assets/6c6048ab-e56c-4ff3-9067-c508047e6c5a)

![MA_UA1061](https://github.com/user-attachments/assets/6d02606f-6ebb-402a-ac0f-d5b327f34392)

![pvalues_UA681](https://github.com/user-attachments/assets/b500d1e7-e995-4ea5-adf2-dd487c545485)

![pvalues_UA1061](https://github.com/user-attachments/assets/93df7ad9-81ff-48c5-afff-779f3929586f)

![fdr_UA681](https://github.com/user-attachments/assets/a5b634f8-1eb4-42d3-ba38-f5a9f3b4751f)

![fdr_UA1061](https://github.com/user-attachments/assets/f6e5f2dd-a74b-47d2-afaf-e4780f590650)

![volcano_UA681](https://github.com/user-attachments/assets/f820e47e-f84e-495b-bff1-342f5ced2fbd)

![volcano_UA1061](https://github.com/user-attachments/assets/535a7bcd-f2fa-48f0-8a48-d674649698c4)

![heatmap](https://github.com/user-attachments/assets/01985909-169f-4d30-aa20-e427dd69bbce)

![scatterplot_log2fc](https://github.com/user-attachments/assets/d86916d2-8053-47bd-848b-023d35964ce0)
