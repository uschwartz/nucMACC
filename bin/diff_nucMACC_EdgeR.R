#!/usr/bin/env Rscript
library(stringr)
library(edgeR)
library(plyranges)
#input <- read.csv("/Volumes/Storage4/THP1_nucMACC/Analysis/Diff_analysis_rep/input_rep.csv")
#path_counts <-"/Volumes/Storage4/THP1_nucMACC/Analysis/Diff_analysis_rep/RUN/02_DIFF_NUCS_READ_COUNTS/sub/diffNucs_readCounts.csv"
#count_summary<-read.delim("/Volumes/Storage4/THP1_nucMACC/Analysis/Diff_analysis_rep/RUN/02_DIFF_NUCS_READ_COUNTS/sub/diffNucs_readCounts.csv.summary")
args   <- commandArgs(TRUE)
input<-read.csv(args[1])
path_counts<-args[2]
count_summary<-read.delim(args[3])
countsmx <- read.table(path_counts, skip = 1, header = T, row.names = 1)
counts<-countsmx[,6:ncol(countsmx)]

coldata <- data.frame(MNase_U = input$MNase_U, Type =input$Type, Rep = input$replicate)
coldata$MNase_U<- as.numeric(coldata$MNase_U)
coldata$Type<-factor(coldata$Type)
coldata$Rep<-factor(coldata$Rep)
rownames(coldata)<-colnames(counts)

DEGL<-DGEList(counts, group = coldata$Type)
DEGL$samples$lib.size <- apply(count_summary[,c(-1)],2,sum)
DEGL<-calcNormFactors(DEGL, "none")
normcounts<-cpm(DEGL, normalized.lib.sizes = T)
design <- model.matrix(~MNase_U + Type + Rep + MNase_U:Type,data = coldata)

y <- estimateGLMCommonDisp(DEGL, design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)

lrt<-glmLRT(fit)
resOrdered<-topTags(lrt, n = Inf)

fdr_cutOFF <- max(which(resOrdered$table$FDR <= 0.05))

positions <- data.frame(chr = str_split_i(rownames(resOrdered$table),"\\.",1), start = as.numeric(str_split_i(rownames(resOrdered$table),"\\.",2)), end = as.numeric(str_split_i(rownames(resOrdered$table),"\\.",3)))
result_table <- cbind(positions, pvalue = resOrdered$table$PValue, padj = resOrdered$table$FDR)

write.table(result_table, "Diff_nucMACC_result_table.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
gr<-as_granges(data.frame(seqnames = result_table$chr, start = result_table$start, end = result_table$end, score = result_table$pvalue))


if(isEmpty(which(resOrdered$table$FDR <= 0.05))==T){
  fdr_cutOFF<-0
  result_table[1,]<-0
}

write.table(result_table[1:fdr_cutOFF,], "Diff_nucMACC_fdrFilt_result_table.tsv", sep = "\t", row.names = F, quote = F, col.names = T)

write_bed(gr[1:fdr_cutOFF,], "Positions_Diff_nucMACC_fdrFilt.bed")
write_bed(gr,"Pvalue_scored_positions.bed")
