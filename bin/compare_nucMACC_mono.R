#!/usr/bin/env Rscript

library(plyranges)
# args[1] = "input.csv"
# args[2] = "Positions_Diff_nucMACC_fdrFilt.bed"
# args[3] = "Diff_nucMACC_fdrFilt_result_table.tsv"
args   <- commandArgs(TRUE)

input <- read.csv(args[1])
Type <- unique(input$Type)
path_analysis <- unique(input$path_analysis)
#path_analysis<-gsub("media/linuxmac","Volumes",path_analysis)
#check input
if(length(Type)!= length(path_analysis)){
  print("There must be one analysis path for each Type!")
}
#load data
diff_mac_nucs <- read_bed(args[2])

nucmacc_1 <- read_bed_graph(paste0(path_analysis[1],"/RUN/05_nucMACC/nucMACC_scores.bedgraph"))
nucmacc_1 <- nucmacc_1[findOverlaps(diff_mac_nucs, nucmacc_1,minoverlap=100)@to]

nucmacc_2 <- read_bed_graph(paste0(path_analysis[2],"/RUN/05_nucMACC/nucMACC_scores.bedgraph"))
nucmacc_2 <- nucmacc_2[findOverlaps(diff_mac_nucs, nucmacc_2,minoverlap=100)@to]

ovrlps <- findOverlapPairs(nucmacc_1, nucmacc_2, minoverlap=100)

acc_down <- ovrlps[which(ovrlps@first$score > ovrlps@second$score)]
diff_down <- acc_down@first$score-acc_down@second$score
acc_down <- acc_down[order(diff_down, decreasing = T)]

acc_up <- ovrlps[which(ovrlps@first$score < ovrlps@second$score)]
diff_up <- acc_up@first$score-acc_up@second$score
acc_up <- acc_up[order(diff_up)]

score(diff_mac_nucs[findOverlaps(diff_mac_nucs,acc_down, minoverlap = 100)@from]) <- diff_down[findOverlaps(diff_mac_nucs,acc_down, minoverlap = 100)@to]
score(diff_mac_nucs[findOverlaps(diff_mac_nucs,acc_up, minoverlap = 100)@from]) <- diff_up[findOverlaps(diff_mac_nucs,acc_up, minoverlap = 100)@to]


if(length(which(score(diff_mac_nucs) == 0)) > 0){
  score(diff_mac_nucs[which(score(diff_mac_nucs) == 0)])<-NA
}

mcols(diff_mac_nucs)$name<-paste0("diff_nuc", 1:length(diff_mac_nucs))

#initiate metadata columns
mcols(diff_mac_nucs)$V1<-rep(NA,length(diff_mac_nucs))
mcols(diff_mac_nucs)$V1.S<-rep(NA,length(diff_mac_nucs))

mcols(diff_mac_nucs)$V2<-rep(NA,length(diff_mac_nucs))
mcols(diff_mac_nucs)$V2.S<-rep(NA,length(diff_mac_nucs))
#assign nucmacc scores
mcols(diff_mac_nucs)$V1.S[findOverlaps(diff_mac_nucs, nucmacc_1, minoverlap = 100)@from]<-score(nucmacc_1[findOverlaps(diff_mac_nucs, nucmacc_1, minoverlap = 100)@to])
mcols(diff_mac_nucs)$V2.S[findOverlaps(diff_mac_nucs, nucmacc_2, minoverlap = 100)@from]<-score(nucmacc_2[findOverlaps(diff_mac_nucs, nucmacc_2, minoverlap = 100)@to])
#add classifications
mcols(diff_mac_nucs)$V1[findOverlaps(diff_mac_nucs, nucmacc_1, minoverlap = 100)@from]<-"normal"
mcols(diff_mac_nucs)$V2[findOverlaps(diff_mac_nucs, nucmacc_2, minoverlap = 100)@from]<-"normal"
rm(nucmacc_1, nucmacc_2)

#load extremes
hyper_1 <- read_bed(paste0(path_analysis[1], "/RUN/05_nucMACC/hyperAcc_monoNucs.bed"))
hyper_2 <- read_bed(paste0(path_analysis[2], "/RUN/05_nucMACC/hyperAcc_monoNucs.bed"))
hypo_1 <- read_bed(paste0(path_analysis[1], "/RUN/05_nucMACC/hypoAcc_monoNucs.bed"))
hypo_2 <- read_bed(paste0(path_analysis[2], "/RUN/05_nucMACC/hypoAcc_monoNucs.bed"))

mcols(diff_mac_nucs)$V1[findOverlaps(diff_mac_nucs, hyper_1, minoverlap = 100)@from]<-"hyper-accessible"
mcols(diff_mac_nucs)$V2[findOverlaps(diff_mac_nucs, hyper_2, minoverlap = 100)@from]<-"hyper-accessible"

mcols(diff_mac_nucs)$V1[findOverlaps(diff_mac_nucs, hypo_1, minoverlap = 100)@from]<-"hypo-accessible"
mcols(diff_mac_nucs)$V2[findOverlaps(diff_mac_nucs, hypo_2, minoverlap = 100)@from]<-"hypo-accessible"

#rename columns
names(diff_mac_nucs@elementMetadata@listData)[3:6] <- c(paste0("Class_",Type[1]),
                                                        paste0("Score_", Type[1]),
                                                        paste0("Class_",Type[2]),
                                                        paste0("Score_", Type[2]))
df <- as.data.frame(diff_mac_nucs)
diff_analysis <- read.table(args[3], header=T)
df <- cbind(df, diff_analysis[,4:5])

write.table(df,"mono_nucMACC_comparison_results.tsv" ,quote = F, sep = "\t", row.names = F)

#get nucs that are more or less accessible, assuming hyer-accessible nucs called in only one condition loose accessibility and hypoaccessibles gain accessibility
T1_higher <- df[unique(c(which(is.na(df[,8])==F & ((df[,9] > df[,11]) | is.na(df[,11]) == T & df[,8] == "hyper-accessible")),
                        which(df[,10]=="hypo-accessible" & is.na(df[,9]==T)) ,which(is.na(df[,10])==F & df[,11] < df[,9]))),]

T1_lower <- df[unique(c(which(is.na(df[,8])==F & df[,9] < df[,11]),
                        which(is.na(df[,10])==F & ((df[,11] > df[,9]) | is.na(df[,9]) == T & df[,10] == "hyper-accessible")),
                        which(df[,8]=="hypo-accessible" & is.na(df[,11]==T)),which(is.na(df[,8])==F & df[,11] > df[,9]))),]

write.table(T1_higher[,1:3],paste0("Gain_in_acc_from_",Type[2],"_to_",Type[1],".bed"), quote = F, sep = "\t", row.names = F, col.names = F)
write.table(T1_lower[,1:3],paste0("Loss_of_acc_from_",Type[2],"_to_",Type[1],".bed"), quote = F, sep = "\t", row.names = F, col.names = F)
