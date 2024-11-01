# RNA-seq
Analysis of RNA-seq data of cassava (Manihot esculenta) plants infected with Xanthomonas phaseoli pv. manihotis from Zárate-Chaves et al. (2021).
This analysis is based on the tutorial by the Cebola Lab (https://github.com/CebolaLab/RNA-seq).

## Pre-alignment QC

## Post-alignment QC

## Visualization 

GC-bias.
```ggg```

```
ggg
```
deeptools computeGCBias -b <sample>-sorted.bam --effectiveGenomeSize 3099922541 -g GCA_000001405.15_GRCh38_no_alt_analysis_set.2bit -l 100 --GCbiasFrequenciesFile <sample>.freq.txt  --biasPlot <sample>.biasPlot.pdf

![SRR13313987 biasPlot](https://github.com/user-attachments/assets/851f67de-4733-4109-957e-24003bce5bc8)
