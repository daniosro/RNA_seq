rm(list = ls())
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("topGO")
library(topGO)
library(tibble)
library(ggplot2)

#Here, we use topGO to determine to assign gene ontologies (GO) to the differentially 
#expressed genes (DEGs)

#Set working directory
setwd('/scratch1/dosorior/RNA_seq')
#Read the table with the GO terms
goids <- read.table('GOid_Me.csv', sep = ',', stringsAsFactors = FALSE)
#Make the first row the column names
colnames(goids) <- goids[1,]
goids <- goids[-1, ] 
head(goids)

#Read the table of DEGs for UA1061
deg_UA1061 <- read.table(file='lfc_UA1061.txt')
deg_UA1061 <- tibble::rownames_to_column(deg_UA1061, "Genes")
head(deg_UA1061)
#Read the table of DEGs for UA681
deg_UA681 <- read.table(file='lfc_UA681.txt')
deg_UA681 <- tibble::rownames_to_column(deg_UA681, "Genes")
head(deg_UA681)

#Merge the differentially expressed genes and the GO terms tables
go_deg_UA1061 <- merge(x=goids, y=deg_UA1061, by = "Genes")
head(go_deg_UA1061)
go_deg_UA681 <- merge(x=goids, y=deg_UA681, by = "Genes")
head(go_deg_UA681)

#Substitute ; for spaces in the GO_IDs column
go_deg_UA1061$GO_IDs <- gsub(';'," ",go_deg_UA1061$GO_IDs)
go_deg_UA681$GO_IDs <- gsub(';'," ",go_deg_UA681$GO_IDs)

write.table(go_deg_UA1061, "go_deg_UA1061.csv",sep="\t", quote = FALSE, col.names = F, row.names = F)
write.table(go_deg_UA681, "go_deg_UA681.csv",sep="\t", quote = FALSE, col.names = F, row.names = F)

GOesByID_UA1061 <- readMappings(file = 'go_deg_UA1061.csv') # list of all GO terms in the transcriptome
GOesByID_UA681 <- readMappings(file = 'go_deg_UA681.csv') # list of all GO terms in the transcriptome

head(GOesByID_UA1061)
head(GOesByID_UA681)

#Make a list of names of the DEGs
UA1061_genes <- names(GOesByID_UA1061)
UA681_genes <- names(GOesByID_UA681)

#Reat table with the protein names
me_proteins <- read.table(file='me_proteins.csv', sep=';', quote="", fill = TRUE)
head(me_proteins)
#Make the first row the column names
colnames(me_proteins) <- me_proteins[1,]
me_proteins <- me_proteins[-1, ] 
head(me_proteins)

#Make a matrix for gene presence/absence between the DEGs and the table assigning genes to protein names
#For UA1061
compared_genes_UA1061 <- factor(as.integer(UA1061_genes %in% me_proteins$Genes))
names(compared_genes_UA1061) <- UA1061_genes
head(compared_genes_UA1061)
#For UA681
compared_genes_UA681 <- factor(as.integer(UA681_genes %in% me_proteins$Genes))
names(compared_genes_UA681) <- UA681_genes
head(compared_genes_UA681)

#topGO categorization for UA1061
BP_myGOdata_UA1061 <- new("topGOdata", description="My project", ontology="BP", allGenes=compared_genes_UA1061, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.gene2GO, gene2GO = GOesByID_UA1061)
CC_myGOdata_UA1061 <- new("topGOdata", description="My project", ontology="CC", allGenes=compared_genes_UA1061, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.gene2GO, gene2GO = GOesByID_UA1061)
MF_myGOdata_UA1061 <- new("topGOdata", description="My project", ontology="MF", allGenes=compared_genes_UA1061, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.gene2GO, gene2GO = GOesByID_UA1061)

# BP_resultFisher <- runTest(BP_myGOdata, algorithm = "classic", statistic = "fisher")
# CC_resultFisher <- runTest(CC_myGOdata, algorithm = "classic", statistic = "fisher")
# MF_resultFisher <- runTest(MF_myGOdata, algorithm = "classic", statistic = "fisher")
?runTest
# An enrichment test was done with the algorithm “weight01” and Kolmogorov-Smirnov test to obtain a list of the top ten GO enriched terms in BP, MF and CC.
BP_resultFisher_UA1061 <- runTest(BP_myGOdata_UA1061, algorithm = "weight01", statistic = "fisher")
CC_resultFisher_UA1061 <- runTest(CC_myGOdata_UA1061, algorithm = "weight01", statistic = "fisher")
MF_resultFisher_UA1061 <- runTest(MF_myGOdata_UA1061, algorithm = "weight01", statistic = "fisher")

# Create and print table with enrichment result
BP_Sanarac_R_DS_SS_UA1061 <- GenTable(BP_myGOdata_UA1061, classicFisher = BP_resultFisher_UA1061, topNodes = 10)
CC_Sanarac_R_DS_SS_UA1061 <- GenTable(CC_myGOdata_UA1061, classicFisher = CC_resultFisher_UA1061, topNodes = 10)
MF_Sanarac_R_DS_SS_UA1061 <- GenTable(MF_myGOdata_UA1061, classicFisher = MF_resultFisher_UA1061, topNodes = 10)
BP_Sanarac_R_DS_SS_UA1061$GO <- 'BP'
CC_Sanarac_R_DS_SS_UA1061$GO <- 'CC'
MF_Sanarac_R_DS_SS_UA1061$GO <- 'MF'
join_fisher_UA1061 <- rbind(BP_Sanarac_R_DS_SS_UA1061, CC_Sanarac_R_DS_SS_UA1061, MF_Sanarac_R_DS_SS_UA1061)

join_fisher_UA1061$treatment <- 'Wilson-DS-Circ'

join_fisher1_UA1061 <- join_fisher_UA1061

join_fisher1_UA1061 <- rbind(join_fisher1_UA1061, join_fisher_UA1061)


join_fisher1_UA1061$classicFisher <- as.numeric(join_fisher1_UA1061$classicFisher)
join_fisher2_UA1061 <- join_fisher1_UA1061
join_fisher2_UA1061$log_p_value <- (-log10(join_fisher1_UA1061$classicFisher))
colnames(join_fisher2_UA1061)

ggplot(data=join_fisher2_UA1061, aes(x= log_p_value, y= reorder(Term, -(log_p_value)), fill=GO)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8) + 
  theme_classic(base_size = 16) + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.2)) + 
  labs(y ="", x="-log10 p-value", fill = 'Tissue') + 
  scale_fill_brewer(palette='Set2') + 
  facet_grid(GO ~ ., scales = "free_y", space = "free_y", margins = F, shrink = T) + 
  theme(legend.position = "none")

write.table(join_fisher2_UA1061, 'DEG_GO_UA1061.txt', sep="\t", quote = FALSE, col.names = T, row.names = F)

#topGO categorization for UA681
BP_myGOdata_UA681 <- new("topGOdata", description="My project", ontology="BP", allGenes=compared_genes_UA681, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.gene2GO, gene2GO = GOesByID_UA681)
CC_myGOdata_UA681 <- new("topGOdata", description="My project", ontology="CC", allGenes=compared_genes_UA681, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.gene2GO, gene2GO = GOesByID_UA681)
MF_myGOdata_UA681 <- new("topGOdata", description="My project", ontology="MF", allGenes=compared_genes_UA681, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.gene2GO, gene2GO = GOesByID_UA681)

# BP_resultFisher <- runTest(BP_myGOdata, algorithm = "classic", statistic = "fisher")
# CC_resultFisher <- runTest(CC_myGOdata, algorithm = "classic", statistic = "fisher")
# MF_resultFisher <- runTest(MF_myGOdata, algorithm = "classic", statistic = "fisher")
?runTest
# An enrichment test was done with the algorithm “weight01” and Kolmogorov-Smirnov test to obtain a list of the top ten GO enriched terms in BP, MF and CC.
BP_resultFisher_UA681 <- runTest(BP_myGOdata_UA681, algorithm = "weight01", statistic = "fisher")
CC_resultFisher_UA681 <- runTest(CC_myGOdata_UA681, algorithm = "weight01", statistic = "fisher")
MF_resultFisher_UA681 <- runTest(MF_myGOdata_UA681, algorithm = "weight01", statistic = "fisher")

# Create and print table with enrichment result
BP_Sanarac_R_DS_SS_UA681 <- GenTable(BP_myGOdata_UA681, classicFisher = BP_resultFisher_UA681, topNodes = 10)
CC_Sanarac_R_DS_SS_UA681 <- GenTable(CC_myGOdata_UA681, classicFisher = CC_resultFisher_UA681, topNodes = 10)
MF_Sanarac_R_DS_SS_UA681 <- GenTable(MF_myGOdata_UA681, classicFisher = MF_resultFisher_UA681, topNodes = 10)
BP_Sanarac_R_DS_SS_UA681$GO <- 'BP'
CC_Sanarac_R_DS_SS_UA681$GO <- 'CC'
MF_Sanarac_R_DS_SS_UA681$GO <- 'MF'
join_fisher_UA681 <- rbind(BP_Sanarac_R_DS_SS_UA681, CC_Sanarac_R_DS_SS_UA681, MF_Sanarac_R_DS_SS_UA681)

join_fisher_UA681$treatment <- 'Wilson-DS-Circ'

join_fisher1_UA681 <- join_fisher_UA681

join_fisher1_UA681 <- rbind(join_fisher1_UA681, join_fisher_UA681)


join_fisher1_UA681$classicFisher <- as.numeric(join_fisher1_UA681$classicFisher)
join_fisher2_UA681 <- join_fisher1_UA681
join_fisher2_UA681$log_p_value <- (-log10(join_fisher1_UA681$classicFisher))
colnames(join_fisher2_UA681)

pdf("plots/GO_UA681.pdf")
ggplot(data=join_fisher2_UA681, aes(x= log_p_value, y= reorder(Term, -(log_p_value)), fill=GO)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8) + 
  theme_classic(base_size = 16) + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.2)) + 
  labs(y ="", x="-log10 p-value", fill = 'Tissue') + 
  scale_fill_brewer(palette='Set2') + 
  facet_grid(GO ~ ., scales = "free_y", space = "free_y", margins = F, shrink = T) + 
  theme(legend.position = "none")
dev.off()

write.table(join_fisher2_UA681, 'DEG_GO_UA681.txt', sep="\t", quote = FALSE, col.names = T, row.names = F)
