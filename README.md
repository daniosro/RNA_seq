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
```

## Pre-alignment Quality Control (QC)

For the quality control assessment, the first step is to activate the previously created environment with the installed packages:

```
#Activate environment
conda activate RNA-seq
```

The QC assessment involves the following steps, described in detail in the Cebola Lab tutorial:

1. Generating a QC report for each sample
2.  Extracting the total number of reads from the QC report
3.  Trimming
4.  Repeating the QC report with fastQC to assess sequence quality after trimming

## Alignment to the reference genome

We downloaded the reference genome for *Manihot esculenta* from https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003957885.1/ and indexed it with STAR.

After indexing, we used Slurm to carry out the alignment from an HPC cluster using STAR. 

## Post-alignment QC

These steps are performed to determine the quality of the genome alignment:

1. Generating QC reports using ```qualimap``` and ```samtools```
2. Combining the QC reports with ```qualimap multi-bamqc```

Below we include two example plots from the output of ```multi-bamqc```, where we can see the genome fraction coverage of each sample and the distribution of sizes of the inserts.

![Genome Fraction Coverage](https://github.com/user-attachments/assets/4f5d8457-3527-4b1a-bf4f-bc613ae5fcdd)

![Insert Size Histogram](https://github.com/user-attachments/assets/a15d8fe3-b7a7-40fa-9cec-83adf0d13572)

## Visualization 

We computed the guanine-citosine (GC) bias of the sequences, which may happen during PCR amplification due to the preferential amplification of specific DNA fragments. We used ```samtools computeGCBias```. We observed a bimodal distribution of the GC bias in all of our samples. Representative examples for 1 replicate of each treatment are included below:

Control

![SRR13313987 biasPlot](https://github.com/user-attachments/assets/851f67de-4733-4109-957e-24003bce5bc8)

UA681

![SRR13313993 biasPlot](https://github.com/user-attachments/assets/04aad013-3540-4b65-8083-f45b8e329867)

UA1061

![SRR13313990 biasPlot](https://github.com/user-attachments/assets/47cfec77-48c1-45e8-8e84-05ee68f221b2)

The bias was corrected with ```correctGCBias```, which removes reads from regions with greater than expected coverage (GC-rich regions) and adds reads from regions with less-than-expected coverage (AT-rich regions).

## Quantification

```Salmon``` was used to generate a transcriptome from the genome files and a matrix of gene counts from the star alignment to the transcriptome. 

Salmon is here used with the variational Bayesian expectation minimisation (VSEM) algorithm for quantification. Quanitifcation is described in the 2020 paper by Deschamps-Francoeur et al., which describes the handling of multi-mapped reads in RNA-seq data. Here, Salmon is run without any normalisation, on each technical replicate; samples are combined and normalised in the next steps.

## Differential expression

All the following code is run in R, and is included in the script RNA-seq_analysis.R

The differential expression analysis contains the following steps:

* Import count data
* Import data to DEseq2
* Determining the differential gene expression
* Generating QC plots before moving on to functional analysis

### QC plots

> PCA
 
We generated a principal component analysis (PCA) plot to determine if the biological replicates clustered together:

![PCA](https://github.com/user-attachments/assets/6fc0c530-924b-452c-a8d7-0a040d751789)

As we can see from the PCA plot, sample 1, which corresponds to the control, and sample 9, which corresponds to a biological replicate of the inoculation with UA681, are outliers, so these samples were excluded from the subsequent analyses. We did not remove any batch effects since all samples corresponded to the same batch.

> MA plots

Next, we generated MA plots (scatter plot of the log fold-change between a sample against the average gene expression (mean of normalised counts)) for each of the treatments:

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
