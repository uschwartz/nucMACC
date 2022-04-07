#!/usr/bin/env Rscript

#### PACKAGE LOADING ####
library(RColorBrewer)

######## DATA LOADING ###########


#defining color palette
dark2 <- c(RColorBrewer::brewer.pal(8, "Dark2"),RColorBrewer::brewer.pal(8, "Set1"),RColorBrewer::brewer.pal(8, "Set2"))


#setting work directory
#work.path <- args[1]
#work.path<-"~/Analysis/MNase_Yeast/H2B_IP/work_H2B_IP_rep1/88/"
#setwd(work.path)


files<-list.files(path="." ,pattern = "insert_size_histogram.txt", full.names = TRUE, recursive = TRUE)

#sample names
samples<-strsplit(files, split = "/")
names<-NULL
for(i in 1:length(samples)){
names<-c(names,samples[[i]][2])
}


#import data
data<-lapply(files, read.table)

#normalize data
for(i in 1:length(data)){
  data[[i]]$V3<-data[[i]]$V2/sum(data[[i]]$V2)
}


#find max value for ylim 
for(i in 1:length(data)){
  curr<-max(data[[i]]$V2)
  if(i==1){
    high<-curr
  }
  else{
    if(curr>high){
      high<-curr
    }
  }
}


#find max value for ylim fraction
for(i in 1:length(data)){
  curr_fr<-max(data[[i]]$V3)
  if(i==1){
    high_fr<-curr_fr
  }
  else{
    if(curr_fr>high_fr){
      high_fr<-curr_fr
    }
  }
}



## show InsertSizeHistogram

pdf("InsertSizeHistogram.pdf",
        width = 8, height = 5)
par(mar=c(6.1, 5.1, 5.1, 18.1),xpd=TRUE)
for(i in 1:length(data)){
  if(i==1){
    plot(data[[i]]$V1,data[[i]]$V2, type="l", col=dark2[i], ylim=c(0,high), xlim=c(55,500), ylab="Counts",xlab="Insert size [bp]", lwd=2)
  }
  else{
    lines(data[[i]]$V1,data[[i]]$V2, type="l", col=dark2[i],lwd=2)
  }
}
par(mar=c(6.1, 5.1, 5.1, 18.1),xpd=FALSE)
abline(v=c(140,200), col="black", lty=2)
par(mar=c(6.1, 5.1, 5.1, 18.1),xpd=TRUE)
legend(550,1.1*high,legend=names, col=dark2[1:length(data)],lty = 1, cex = 0.8, bty="n", lwd=2)
dev.off()

pdf("InsertSizeHistogram_Fraction.pdf",
    width = 8, height = 5)
par(mar=c(6.1, 5.1, 5.1, 18.1),xpd=TRUE)
for(i in 1:length(data)){
  if(i==1){
    plot(data[[i]]$V1,data[[i]]$V3, type="l", col=dark2[i], ylim=c(0,high_fr),xlim=c(55,500), ylab="Fraction",xlab="Insert size [bp]",lwd=2)
  }
  else{
    lines(data[[i]]$V1,data[[i]]$V3, type="l", col=dark2[i], lwd=2)
  }
}
par(mar=c(6.1, 5.1, 5.1, 18.1),xpd=FALSE)
abline(v=c(140,200), col="black", lty=2)
par(mar=c(6.1, 5.1, 5.1, 18.1),xpd=TRUE)
legend(550,1.1*high_fr,legend=names, col=dark2[1:length(data)],lty = 1, cex = 1, bty="n", lwd=2)
dev.off()
