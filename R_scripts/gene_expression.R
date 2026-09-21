library(DESeq2)
library(tidyverse)
library(rlist)
library(readr)
library(dplyr)
#Compare Ht vs Wt gene
#read in data

Counts1 <- as.matrix (read.csv("gene_count_matrix.csv", row.names = "gene_id"))

Coldata <- read.csv("ColData.csv", row.names = 1)


dds <- DESeqDataSetFromMatrix(countData = Counts1,
                              colData = Coldata,
                              design= ~ Condition)
dds <- DESeq(dds)
#remove genes/rows with less than 10 transcript count
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Compare HT with WT (drug)
HTvsWT_drug<-results(dds, contrast = c("Condition", "HTdrug","WTdrug")) 
HTvsWT_drug <- na.omit(HTvsWT_drug)



#Filter for , logFold2, pvalue 
 
HTvsWT_drug <- HTvsWT_drug[(HTvsWT_drug$log2FoldChange > 0.58 | HTvsWT_drug$log2FoldChange < -0.58) & HTvsWT_drug$padj <= 0.05, ]
write.csv(HTvsWT_drug, "DEGs_HTvsWT_drug.csv")

#compare HT, WT (veh)
HTvsWT_veh<-results(dds, contrast = c("Condition", "HTveh","WTveh")) 
HTvsWT_veh <- na.omit(HTvsWT_veh)

#Filter for , logFold2, pvalue 


HTvsWT_veh <- HTvsWT_veh[(HTvsWT_veh$log2FoldChange > 0.58 | HTvsWT_veh$log2FoldChange < -0.58) & HTvsWT_veh$padj <= 0.05, ]
write.csv(HTvsWT_drug, "DEGs_HTvsWT_veh.csv")


#for normalized gene counts

normalized_counts_genes <- counts(dds, normalized=TRUE)
write.csv(normalized_counts_genes, "normalised_ADAR_counts_gene.csv")

##############################################

#Transcript expression

#read in data

Counts2 <- as.matrix (read.csv("CountMatrix/transcript_count_matrix.csv", row.names = "transcript_id"))

Coldata <- read.csv("CountMatrix/ColData.csv", row.names = 1)


dds <- DESeqDataSetFromMatrix(countData = Counts2,
                              colData = Coldata,
                              design= ~ Condition)
dds <- DESeq(dds)
#remove genes/rows with less than 10 transcript count
keep <- rowSums(counts(dds,)) >= 10
dds <- dds[keep,]

#Compare Ht with Wt for vehicle
HtvsWt_transcript <-results(dds, contrast = c("Condition", "Ht","Wt")) 
HtvsWt_transcript <- na.omit(HtvsWt_transcript)
#Filter for P_value =< 0.05, logFold2 
HtvsWt_transcript_pvalue<- HtvsWt_transcript[(HtvsWt_transcript$log2FoldChange > 0.58 | HtvsWt_transcript$log2FoldChange < -0.58) & HtvsWt_transcript$pvalue <= 0.05, ]
print(HtvsWt_pvalue)
#filter for p_adj
HtvsWt_transcript_padj <- HtvsWt_transcript[(HtvsWt_transcript$log2FoldChange > 0.58 | HtvsWt_transcript$log2FoldChange < -0.58) & HtvsWt_transcript$padj <= 0.05, ]
print(HtvsWt_padj)

#save results
write.csv(HtvsWt_transcript_pvalue, "HtvsWt_transcript_pvalue.csv")
write.csv(HtvsWt_transcript_padj, "HtvsWt_transcript_padj.csv")

#retrieve normalized counts
normalized_counts_transcripts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts_transcripts, "normalized_counts_transcripts.csv")
#for genes
normalized_counts_genes <- counts(dds, normalized=TRUE)
write.csv(normalized_counts_genes, "normalised_counts_gene.csv")


















