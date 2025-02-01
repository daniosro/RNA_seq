# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

update.packages()
#BiocManager::install("cqn") #optional for cqn normalisation
BiocManager::install("DESeq2")
a# BiocManager::install("tximport")
# BiocManager::install("biomaRt")
BiocManager::install("GenomicFeatures", force = TRUE)
# BiocManager::install("txdbmaker")
# install.packages("Rcpp", type = "source")
# BiocManager::install("apeglm")
# BiocManager::install("pheatmap")

.rs.restartR()

remove.packages(DESeq2)
library(BiocManager)
library(GenomicFeatures)
library(txdbmaker)
library(tximport)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(gplots)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(RColorBrewer)
library(dplyr)

setwd('..')
setwd('..')
setwd('/scratch1/dosorior/RNA_seq/')

#Read in the files with the sample information
samples = read.table('samples3.txt')
# Make a transcript to gene ID file
txdb <- makeTxDbFromGFF("gca_annotation_no_genes.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
write.table(tx2gene,"tx2gene.txt",sep="\t")
#Read the count data with tximport
counts.imported = tximport(files = as.character(samples[,2]), type = 'salmon', tx2gene = tx2gene)

#Create dataframe with sample information for DEseq2
colData <- data.frame(
  condition = c("Control", "Control", "Control", "UA 1061", "UA 1061", "UA 1061", "UA681", "UA681", "UA681"),
  row.names = c("Control_1", "Control_2", "Control_3", "UA1061_1", "UA1061_2", "UA1061_3", "UA681_1", "UA681_2", "UA681_3"),
  stringsAsFactors = FALSE
)
colData$condition = factor(colData$condition)

#Import to DEseq2
counts.DEseq = DESeqDataSetFromTximport(counts.imported, colData = colData, design = ~condition)

dds <- DESeq(counts.DEseq)
resultsNames(dds) #lists the coefficients

plotDispEsts(dds)

#Determine differential gene expression
# Load the apeglm library to apply the apeglm shrinkage method to shrink high log-fold 
#changes with little statistical evidence and account for lowly expressed genes with significant deviation

#List the names of the coefficients and choose your comparison
resultsNames(dds)

#Apply the LFC shrinkage
#For UA681 vs control
res_lfc_UA681 <- lfcShrink(dds, coef=3, type="apeglm")
#For UA1061 vs control
res_lfc_UA1061 <- lfcShrink(dds, coef=2, type="apeglm")

#Quality control
#Principal component analysis
#For a few samples like this case
rld <- rlog(dds, blind = TRUE)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

#Plotting the results with ggplot

z = plotPCA(rld, "condition")
nudge <- position_nudge(y = 2,x=6)
z = z + geom_text(aes(label = name), position = nudge) + theme()
pdf("plots/PCA.pdf")
print(z)
dev.off()

#plotPCA from DEseq2 plots uses the top 500 genes:
data = plotPCA(rld, intgroup = ("condition"), returnData = TRUE)
p <- ggplot(data, aes(x = PC1, y = PC2, color = condition ))
p <- p + geom_point() + theme ()
print(p)

#Alternatively, PCA can be carried out using all genes:
df_out <- as.data.frame(pca$x)
df_out$group <- samples[,3]

#Include the next two lines to add the PC % to the axis labels
percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), paste0(" (", as.character(percentage), "%", ")"), sep="")

p <- ggplot(df_out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point() + theme() + xlab(percentage[1]) + ylab(percentage[2])

print(p)


#Generate MA plot (scatter plot of the log fold-change between 
#two samples against the average gene expression (mean of normalised counts)).

plotMA(res_lfc_UA681, main = 'UA681', cex = 0.5)
plotMA(res_lfc_UA1061, main = 'UA1061', cex = 0.5)

#Since the PCA shows that control 1 is clustering with the samples, and than replicate 9 is clustering separately
#from all the other samples, these 2 will be excluded from further analysis

#Removing columns 1 and 9 from the countdata dataframe
dds_1 <- dds[,-1]
dds_clean <- dds_1[,-8]
dds_clean_final<-DESeq(dds_clean)

resultsNames(dds_clean_final) #lists the coefficients
pdf(file="plots/DispEsts.pdf")
plotDispEsts(dds_clean_final)
dev.off()

#Apply the LFC shrinkage
#For UA681 vs control
res_lfc_UA681 <- lfcShrink(dds_clean_final, coef=3, type="apeglm")
#For UA1061 vs control
res_lfc_UA1061 <- lfcShrink(dds_clean_final, coef=2, type="apeglm")

#Quality control
#Principal component analysis
#For a few samples like this case
rld <- rlog(dds_clean_final, blind = TRUE)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

#Plotting the results with ggplot
z = plotPCA(rld, "condition")
nudge <- position_nudge(y = 2,x=6)
z = z + geom_text(aes(label = name), position = nudge) + theme()
pdf("plots/PCA_filtered.pdf")
print(z)
dev.off()

#Fancier PCA plot with all genes:
samples_filtered <- c("Control", "Control", "UA1061", "UA1061", "UA1061", "UA681", "UA681")
df_out <- as.data.frame(pca$x)
df_out$group <- samples_filtered

#Include the next two lines to add the PC % to the axis labels
percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
percentage <- paste(colnames(df_out), paste0(" (", as.character(percentage), "%", ")"), sep="")

p <- ggplot(df_out, aes(x = PC1, y = PC2, color=group))
p <- p + geom_point() + theme() + xlab(percentage[1]) + ylab(percentage[2])

print(p)


#Generate MA plot (scatter plot of the log fold-change between 
#two samples against the average gene expression (mean of normalised counts)).

pdf("plots/MA_UA681.pdf")
plotMA(res_lfc_UA681, main = 'UA681', cex = 0.5)
dev.off()

pdf("plots/MA_UA1061.pdf")
plotMA(res_lfc_UA1061, main = 'UA1061', cex = 0.5)
dev.off()

#Distribution of p-values and FDRs
#For UA681
#The distribution of p-values
pdf("plots/pvalues_UA681.pdf")
hist(res_lfc_UA681$pvalue, breaks = 50, col = 'grey', main = 'UA681', xlab = 'p-value')
dev.off()
#The false-discovery rate distribution
pdf("plots/fdr_UA681.pdf")
hist(res_lfc_UA681$padj, breaks = 50, col = 'grey', main = 'UA681', xlab = 'FDR')
dev.off()

#For UA1061
#The distribution of p-values
pdf("plots/pvalues_UA1061.pdf")
hist(res_lfc_UA1061$pvalue, breaks = 50, col = 'grey', main = 'UA1061', xlab = 'p-value')
dev.off()
#The false-discovery rate distribution
pdf("plots/fdr_UA1061.pdf")
hist(res_lfc_UA1061$padj, breaks = 50, col = 'grey', main = 'UA1061', xlab = 'FDR')
dev.off()

#Volcano plots
#A volcano plot is a scatterplot which plots the p-value of differential expression 
#against the fold-change. The volcano plot can be designed to highlight datapoints of 
#significant genes, with a p-value and fold-change cut off.

#Allow for more space around the borders of the plot
par(mar = c(5, 4, 4, 4))

#Set your log-fold-change and p-value thresholds
lfc = 2
pval = 0.05

#For UA681
tab = data.frame(logFC = res_lfc_UA681$log2FoldChange, negLogPval = -log10(res_lfc_UA681$padj))#make a data frame with the log2 fold-changes and adjusted p-values

pdf('plots/volcano_UA681.pdf')
p <- plot(tab, pch = 16, cex = 0.4, xlab = expression(log[2]~fold~change),
          ylab = expression(-log[10]~pvalue), main = 'UA681') 

#Genes with a fold-change greater than 2 and p-value<0.05:
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))

#Colour these red
p <- p + points(tab[signGenes, ], pch = 16, cex = 0.5, col = "red")

#Show the cut-off lines
p <- p + abline(h = -log10(pval), col = "green3", lty = 2)
p <- p + abline(v = c(-lfc, lfc), col = "blue", lty = 2)

#Add extra text
p <- p + mtext(paste("FDR =", pval), side = 4, at = -log10(pval), cex = 0.6, line = 0.5, las = 1)

p <- p + mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),
               cex = 0.6, line = 0.5)
print(p)
dev.off()

#For UA1061
tab = data.frame(logFC = res_lfc_UA1061$log2FoldChange, negLogPval = -log10(res_lfc_UA1061$padj))#make a data frame with the log2 fold-changes and adjusted p-values

pdf('plots/volcano_UA1061.pdf')
plot(tab, pch = 16, cex = 0.4, xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue), main = 'UA1061') 

#Genes with a fold-change greater than 2 and p-value<0.05:
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))

#Colour these red
points(tab[signGenes, ], pch = 16, cex = 0.5, col = "red")

#Show the cut-off lines
abline(h = -log10(pval), col = "green3", lty = 2)
abline(v = c(-lfc, lfc), col = "blue", lty = 2)

#Add extra text
mtext(paste("FDR =", pval), side = 4, at = -log10(pval), cex = 0.6, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),
      cex = 0.6, line = 0.5)
print(p)
dev.off()

#Heatmap

#Make a z-score heat map
res_dds <- results(dds_clean_final)
diffgenes <- rownames(res_dds)[ which(res_dds$padj < 0.05) ]
normCounts <- counts(dds_clean_final, normalized = TRUE)
diffcounts <- normCounts[ diffgenes, ]

my_palette <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(diffcounts, 
          labRow = "", dendrogram="row",
          trace = "none", density.info = "none",
          scale = "row", col=my_palette,
          symm=F,symkey=T,symbreaks=T,
          distfun = function(x) as.dist(1 - cor(t(x))))

#For log2 fold change
#For UA681
top_UA681 <- res_lfc_UA681[ (res_lfc_UA681$baseMean > 50) & (abs(res_lfc_UA681$log2FoldChange) > 2),] 
diffgenes_UA681 <- rownames(top_UA681)
pp <- res_lfc_UA681[diffgenes_UA681,]$log2FoldChange
l2_val_UA681 <- as.matrix(res_lfc_UA681[diffgenes_UA681,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val_UA681)<-"logFC"
#For UA1061
top_UA1061 <- res_lfc_UA1061[ (res_lfc_UA1061$baseMean > 50) & (abs(res_lfc_UA1061$log2FoldChange) > 2),] 
diffgenes_UA1061 <- rownames(top_UA1061)
l2_val_UA1061 <- as.matrix(res_lfc_UA1061[diffgenes_UA1061,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val_UA1061)<-"logFC"


#maps values between b/w/r for min and max l2 values
col_logFC_UA681 <- colorRamp2(c(-10,0, 10), c("red", "white", "blue")) 
col_logFC_UA1061 <- colorRamp2(c(min(l2_val_UA1061),0, max(l2_val_UA1061)), c("red", "white", "blue")) 

merged <- merge(l2_val_UA681,l2_val_UA1061,by="row.names",all.x=TRUE)
merged[is.na(merged)] <- 0 
merged[1]

#Function to make a matrix from dataframe
matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}
#Convert log2 dataframe to matrix
log2matrix <- matrix.please(merged)
#Replace names of columns
colnames(log2matrix) <- c('UA681','UA1061')

#Make heatmap
ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))
pdf("plots/heatmap.pdf")
h <- Heatmap(log2matrix, 
             cluster_rows = T, name="Log2FC", 
             cluster_columns = F, col = col_logFC_UA681,
             show_row_names = FALSE,
             row_dend_width = unit(4, "cm")
)
draw(h)
dev.off()


#Build correlation plot
#Extracting pvalues for all the genes
#For UA681
p_val_UA681 <- as.matrix(res_lfc_UA681$pvalue) #getting  pvalue for each gene 
rownames(p_val_UA681) <- rownames(res_lfc_UA681)
colnames(p_val_UA681) <- "pvalue"
#For UA1061
p_val_UA1061 <- as.matrix(res_lfc_UA1061$pvalue) #getting  pvalue for each gene 
rownames(p_val_UA1061) <- rownames(res_lfc_UA1061)
colnames(p_val_UA1061) <- "pvalue"

#Extracting log2FC for all the genes 
#For UA681
log2fc_UA681 <- as.matrix(res_lfc_UA681$log2FoldChange) #getting  pvalue for each gene 
rownames(log2fc_UA681) <- rownames(res_lfc_UA681)
colnames(log2fc_UA681) <- "log2FC"
#For UA1061
log2fc_UA1061 <- as.matrix(res_lfc_UA1061$log2FoldChange) #getting  pvalue for each gene 
rownames(log2fc_UA1061) <- rownames(res_lfc_UA1061)
colnames(log2fc_UA1061) <- "log2FC"

#Extracting BaseMean for all the genes 
#For UA681
baseMean_UA681 <- as.matrix(res_lfc_UA681$baseMean) #getting  pvalue for each gene 
rownames(baseMean_UA681) <- rownames(res_lfc_UA681)
colnames(baseMean_UA681) <- "baseMean"
#For UA1061
baseMean_UA1061 <- as.matrix(res_lfc_UA1061$baseMean) #getting  pvalue for each gene 
rownames(baseMean_UA1061) <- rownames(res_lfc_UA1061)
colnames(baseMean_UA1061) <- "baseMean"


#Merge pvalues from both treatments in a single dataframe
merged_pval <- merge(p_val_UA681, p_val_UA1061,by="row.names",all.x=TRUE)
#Merge log2FC from both treatments in a single dataframe
merged_log2fc <- merge(log2fc_UA681, log2fc_UA1061,by="row.names",all.x=TRUE)
#Merge baseMean from both treatments in a single dataframe
merged_baseMean <- merge(baseMean_UA681, baseMean_UA1061,by="row.names",all.x=TRUE)
#Merge Log2FC and pvalues dataframes
merged_pval_log2fc <- merge(merged_log2fc, merged_pval,by="row.names",all.x=TRUE)
#Remove unnecesary columns from dataframe after merge
#Check columns to remove
colnames(merged_pval_log2fc)
#Remove columns
merged_pval_log2fc_filt <- merged_pval_log2fc %>% select(-c(Row.names, Row.names.y))
#Rename row names column
colnames(merged_pval_log2fc_filt)[which(names(merged_pval_log2fc_filt) == "Row.names.x")] <- "Row.names"

#Merge baseMean, Log2FC and pvalues dataframes
merged_final <- merge(merged_baseMean, merged_pval_log2fc_filt,by="row.names",all.x=TRUE)
#Remove unnecesary columns from dataframe after merge
#Check columns to remove
colnames(merged_final)
#Remove columns
merged_final_filt <- merged_final %>% select(-c(Row.names, Row.names.y))
#Replace NAs
merged_final[is.na(merged_final)] <- 0 

#Replace names of columns
colnames(merged_final_filt) <- c('genes', 'baseMean_UA681', 'baseMean_UA1061', 'logFC_UA681','logFC_UA1061', 'pvalue_UA681', 'pvalue_UA1061')

df_scatterplot <- merged_final_filt %>% mutate(category = case_when(
  ((logFC_UA681 < -2 & logFC_UA1061 > 2) & (pvalue_UA681 < 0.05 | pvalue_UA1061 < 0.05)) |
    ((logFC_UA681 > 2 & logFC_UA1061 < -2) & (pvalue_UA681 < 0.05 | pvalue_UA1061 < 0.05))  ~ "Contrasting",
  logFC_UA681 > 2 & logFC_UA1061 > 2 & pvalue_UA681 < 0.05 & pvalue_UA1061 < 0.05 &
    baseMean_UA681 > 50 & baseMean_UA1061 > 50 ~ "Overexpressed",
  logFC_UA681 < -2 & logFC_UA1061 < -2 & pvalue_UA681 < 0.05 & pvalue_UA1061 < 0.05 &
    baseMean_UA681 > 50 & baseMean_UA1061 > 50 ~ "Repressed",
  TRUE ~ "Unchanged"
))
df_scatterplot 

# Make scatterplot
pdf("plots/scatterplot_log2fc.pdf")
ggplot(df_scatterplot, aes(x=logFC_UA681, y=logFC_UA1061, color=category)) + 
  geom_point() +
  coord_fixed(ratio = 1) +
  scale_x_continuous(name="Log2FC (UA681)", breaks=seq(-10,10,10)) +
  scale_y_continuous(name="Log2FC (UA1061)",breaks=seq(-10,10,10)) +
  scale_color_manual(values=c("red","blue", "grey")) +
  theme_bw() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=14,face="bold"))

dev.off()

#Data exploration

#1. How many genes are differentially expressed?

#For UA681
#increased expression
attach(as.data.frame(res_lfc_UA681))

#The total number of DEGs with an adjusted p-value<0.05
summary(res_lfc_UA681, alpha=0.05)

#The total number of DEGs with an adjusted p-value<0.05 AND absolute fold-change > 2
sum(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >2)

#Decreased expression:
sum(!is.na(padj) & padj < 0.05 & log2FoldChange <0) #any fold-change
sum(!is.na(padj) & padj < 0.05 & log2FoldChange <(-2)) #fold-change greater than 2

#Increased expression:
sum(!is.na(padj) & padj < 0.05 & log2FoldChange >0) #any fold-change
sum(!is.na(padj) & padj < 0.05 & log2FoldChange >2) #fold-change greater than 2


#For UA1061
#increased expression
attach(as.data.frame(res_lfc_UA1061))

#The total number of DEGs with an adjusted p-value<0.05
summary(res_lfc_UA1061, alpha=0.05)

#The total number of DEGs with an adjusted p-value<0.05 AND absolute fold-change > 2
sum(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >2)

#Decreased expression:
sum(!is.na(padj) & padj < 0.05 & log2FoldChange <0) #any fold-change
sum(!is.na(padj) & padj < 0.05 & log2FoldChange <(-2)) #fold-change greater than 2

#Increased expression:
sum(!is.na(padj) & padj < 0.05 & log2FoldChange >0) #any fold-change
sum(!is.na(padj) & padj < 0.05 & log2FoldChange >2) #fold-change greater than 2

# What are the top genes?

#At this stage it may be useful to create a copy of the results with the gene version removed from the gene name, to make it easier for you to search for the gene name etc. 
#The rownames currently appear as 'ENSG00000175197.12, ENSG00000128272.15' etc.
#To change them to 'ENSG00000175197, ENSG00000128272'
# For UA681
lfc_UA681.gene = as.data.frame(res_lfc_UA681)

#Some gene names are repeated if they are in the PAR region of the Y chromosome. Since dataframes cannot have duplicate row names, we will leave these gene names as they are and rename the rest.
whichgenes = which(!grepl('PAR', rownames(lfc_UA681.gene)))
rownames(lfc_UA681.gene)[whichgenes] = unlist(lapply(strsplit(rownames(lfc_UA681.gene)[whichgenes], '\\.'), '[[',1))

#remove v8 from rownames
rownames(lfc_UA681.gene) = substr(rownames(lfc_UA681.gene),1,15)
write.table(lfc_UA681.gene, file='lfc_UA681.txt', quote=FALSE)
head(lfc_UA681.gene)
colnames(lfc_UA681.gene)
a <- subset(lfc_UA681.gene, padj < 0.05)
head(a)

# #Decreased expression
# ql53.DEGs.down <- groups12.table$FDR < 0.05 & groups12.table$logFC<0
# names(ql53.DEGs.down) <- rownames(groups12.table)
# pwf.dn <- nullp(ql53.DEGs.up, "hg19", "ensGene")
# go.results.dn <- goseq(pwf.dn, "hg19", "ensGene")
# 
# #Increased expression
# ql53.DEGs.up <- groups12.table$FDR < 0.05 & groups12.table$logFC>0
# names(ql53.DEGs.up) <- rownames(groups12.table)
# pwf.up <- nullp(ql53.DEGs.down, "hg19","ensGene")
# go.results.up <- goseq(pwf.up, "hg19","ensGene")

#subset the significant genes
lfc_UA681.sig = subset(lfc_UA681.gene, padj < 0.05 & !is.na(padj))#subset the significant genes
head(lfc_UA681.sig)
#Export table with significant genes
write.table(lfc_UA681.sig, file='lfc_UA681_sig.txt', quote=FALSE)

#View the top 10 genes with the most significant (adjusted) p-values
head(lfc_UA681.sig, n = 10)
dim(lfc_UA681.sig)
#The largest fold-changes with a significant p-value
lfc_UA681.sig[order(abs(lfc_UA681.sig$log2FoldChange), decreasing = TRUE),][1:10,] #add the [1:10,] to see the top 10 rows
# For UA1061
lfc_UA1061.gene = as.data.frame(res_lfc_UA1061)
rownames(lfc_UA1061.gene)
#Some gene names are repeated if they are in the PAR region of the Y chromosome. Since dataframes cannot have duplicate row names, we will leave these gene names as they are and rename the rest.
whichgenes = which(!grepl('PAR', rownames(lfc_UA1061.gene)))
rownames(lfc_UA1061.gene)[whichgenes] = unlist(lapply(strsplit(rownames(lfc_UA1061.gene)[whichgenes], '\\.'), '[[',1))

#remove v8 from rownames
rownames(lfc_UA1061.gene) = substr(rownames(lfc_UA1061.gene),1,15)
write.table(lfc_UA1061.gene, file='lfc_UA1061.txt', quote=FALSE)

#subset the significant genes
lfc_UA1061.sig = subset(lfc_UA1061.gene, padj < 0.05 & !is.na(padj))#subset the significant genes
dim(lfc_UA1061.sig)
#Export table with significant genes
write.table(lfc_UA1061.sig, file='lfc_UA1061_sig.txt', quote=FALSE)

#View the top 10 genes with the most significant (adjusted) p-values
head(lfc_UA1061.sig, n = 10)

#The largest fold-changes with a significant p-value
lfc_UA1061.sig[order(abs(lfc_UA1061.sig$log2FoldChange), decreasing = TRUE),][1:10,] #add the [1:10,] to see the top 10 rows

#Plot the expression for the top genes
#Select your chosen gene 
tmp = plotCounts(dds_clean_final, gene = grep('MANES_01G000500v8', names(dds_clean_final), value = TRUE), intgroup = "condition", pch = 18, main = '??? expression', returnData = TRUE)

theme <- theme(panel.background = element_blank(), panel.border = element_rect(fill = NA),
               plot.title = element_text(hjust = 0.5))

p <- ggplot(tmp, aes(x = condition, y = count)) + geom_boxplot() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.6) + ggtitle('??? expression') + theme ()

print(p)

