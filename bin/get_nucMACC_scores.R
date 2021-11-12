#!/usr/bin/env Rscript

#### PACKAGE LOADING ####
library(LSD)

######## DATA LOADING ###########
args   <- commandArgs(TRUE)

#read sample input file
input<-read.csv(file=args[1])

## read count table
readCounts <- read.delim(file = args[2],header = T)


#get total counts statistics
featureCounts.mono<-read.delim(file=args[3])
colnames(featureCounts.mono)<-gsub("_mono.bam", "", colnames(featureCounts.mono))

dir.create("Figures", showWarnings = FALSE)

#################################################
################### nucMACC #####################
#################################################
### DATA PREPARATION

#get count table
count.table<-readCounts[,c(6:(ncol(readCounts)-2))]
rownames(count.table)<-readCounts$nucID
colnames(count.table)<-gsub("_mono.bam", "", colnames(count.table))

### READ filtering based on raw reads
#raw reads threshold
raw.flt<-30

#filter
idx.raw<-apply(count.table,1,sum)>raw.flt
realNucs <- count.table[idx.raw,]

### READ COUNT NORMALIZATION
#get normalization factors for CPMs
normFactor<-apply(featureCounts.mono[,c(-1)],2,sum)/1e6

#Counts per million (normalization based on library size)
normCounts<-t(t(realNucs)/normFactor[match(colnames(realNucs),
                                               names(normFactor))]) 



#get corresponding MNase conc
mx<-match(colnames(normCounts), input$Sample_Name)
mnase_conc<-input[mx,"MNase_U"]

## get pseudocount
pseudocount<-median(apply(normCounts,2,median))

#The slope of linear regression
makeRegr <- function(x){
    titration.points <- log2(mnase_conc)
    x.norm<-log2(x+pseudocount)
    fit <- lm(x.norm~titration.points)
    slope <- fit$coefficients[2]*-1
    R2<-summary(fit)$r.squared
    return(data.frame(slope,R2))
}

regr.list <- apply(normCounts, 1, makeRegr)    #Actual calculation of linear regression
regr.results <-do.call(rbind, regr.list)

#####################################################
############## nucMACC GC correction  ###############
#####################################################

#extract nucMACC
slope <- regr.results$slope
names(slope)<-rownames(regr.results)

#add GC
mx.gc<-match(names(slope), as.character(readCounts$nucID))
gc <- readCounts$GC_cont[mx.gc]
names(gc)<-readCounts$nucID[mx.gc]

###########################################
########## LOWESS GC CORRECTION ###########
###########################################

## filter nucs exhibiting extrem  GC content
# quantify the number of nucs within 5% GC content intervals
quant<-hist(gc, breaks=seq(0,1,by = 0.05) ,  plot = F)


# at least 200 nucleosomes required per 5% step
#thresholds
max.gc<-quant$breaks[max(which(quant$counts > 200)+1)]
min.gc<-quant$breaks[min(which(quant$counts > 200))]

#filter based on GC
gc.filt<-gc[which(gc > min.gc & gc < max.gc)]
slope.filt<-slope[which(gc > min.gc & gc < max.gc)]


#calculate loess fit
loess.slope <- loess(y ~ x, span=0.1,
                     data.frame(x=gc.filt, y=slope.filt),
                     control=loess.control(surface = "interpolate",
                                           statistics="none",cell=0.1))


#predict norm. factor for the slope
predict.slope <- predict(loess.slope, data.frame(x=gc.filt))
predict.slope[is.na(predict.slope)] <- 0


#### PLOT CORRELATION BEFORE THE GC NORMALIZATION with linar fit slope
ord <- order(gc.filt)


png("Figures/before_gc_norm_finalPlot.png", width = 1280, height = 1280, res =300 )
    heatscatter(x=gc.filt[ord],y=slope.filt[ord], cor=T, cexplot=0.5, ylab="nucMACC scores",
                xlab="GC content")
    abline(h=0,col="black",lwd=3,lty=2)
    lines(gc.filt[ord],predict.slope[ord], col="#5CB85B",lwd=3)
dev.off()

#### PLOT OF CORRELATION AFTER THE GC NORMALIZATION
png("Figures/after_gc_norm_finalPlot.png", width = 1280, height = 1280, res =300 )
    heatscatter(x=gc.filt[ord], y=(slope.filt[ord])-(predict.slope[ord]-median(slope.filt)),
                cor=T,cexplot=0.5, ylab="nucMACC scores", xlab="GC content")
    abline(h=0,col="black",lwd=3,lty=2)
dev.off()



# GC normalize nucMACC scores
nucMACC_scores <- slope.filt-(predict.slope-median(slope.filt))

########################################################
########## get hypo- and hyper-accessible Nucs #########
########################################################

# prepare for calling
# normalize nucMACC scores to maximum 1
nucMACC_sort <- sort(nucMACC_scores)
nucMACC_norm <- (nucMACC_sort/(max(nucMACC_sort)-min(nucMACC_sort)))
#normalize ranks to maximum 1
nucMACC_rank <- rank(nucMACC_sort, ties.method = "first")/length(nucMACC_sort)

### get LOESS smoothing function
#define to take 1000 nucleosomes as window for loess
span.loess<-1000/length(nucMACC_scores)
lm <- loess(nucMACC_norm~nucMACC_rank,span=span.loess)

#Curve fitting
loess.line <- predict(lm)

#get first derivative to deduce the slope of loess function
loess.f1 <- diff(loess.line)/diff(nucMACC_rank)

#set the slope cutOff to 1
cutOff1 <- min(which(loess.f1 < 1))
cutOff2 <- max(which(loess.f1 < 1))

## show loess derivative and cutoff 1
png("Figures/loess_first_order_derivative.png",
        width = 1280, height = 1280, res =300 )
    plot(loess.f1, xlab="nucMACC rank", ylab="loess first order derivative",
         type="l")
    abline(h = 1, col="#D43F39")
dev.off()

#final selection
png("Figures/nucMACC_selection.png",
    width = 1280, height = 1280, res =300 )
    plot(nucMACC_rank,nucMACC_norm, ylab="normalized nucMACC score",
         xlab = "normalized rank", bty="n", pch=15)
    abline(v=nucMACC_rank[cutOff1], col="#984EA3", lty=2)
    abline(v=nucMACC_rank[cutOff2], col="#1B9E77", lty=2)
    points(nucMACC_rank[1:cutOff1], nucMACC_norm[1:cutOff1],pch=15, col="#984EA3")
    points(nucMACC_rank[cutOff2:length(nucMACC_rank)],
           nucMACC_norm[cutOff2:length(nucMACC_rank)],pch=15, col="#1B9E77")
    lines(nucMACC_rank,loess.line, col="red")
dev.off()

### select the nucs with nucMACC scores deviating from the mean
nucMACC_low <- nucMACC_scores < nucMACC_sort[cutOff1]
nucMACC_high <- nucMACC_scores > nucMACC_sort[cutOff2]


### generate data table
df<-data.frame(nucID=names(nucMACC_scores), nucMACC=nucMACC_scores)

df$category<-"normal"
df$category[nucMACC_high]<-"hyper-accessible"
df$category[nucMACC_low]<-"hypo-accessible"

##concatenate tables
nucStats<-cbind(readCounts[match(df$nucID,readCounts$nucID),
                           c("Chr", "Start","End", "Strand","GC_cont")],df,
                regr.results[match(df$nucID,rownames(regr.results)),])

write.table(nucStats, file="nucMACC_result_table.tsv",
            row.names = FALSE, sep="\t", quote=FALSE,col.names = T)

#export nucMACC as bedgraph
bedgraph <- nucStats[,c("Chr","Start","End","nucMACC")]
write.table(bedgraph , file="nucMACC_scores.bedgraph",
            row.names = FALSE, sep="\t", quote=FALSE, col.names = FALSE)


#export bed file and define score nucMACC
#hypo accessible
hypoAcc<-subset(nucStats, category=="hypo-accessible")
hypoAcc_bed <- hypoAcc[,c("Chr","Start","End","nucID","nucMACC","Strand")]
write.table(hypoAcc_bed , file="hypoAcc_monoNucs.bed",
            row.names = FALSE, sep="\t", quote=FALSE, col.names = FALSE)

#hyper accessible
hyperAcc<-subset(nucStats, category=="hyper-accessible")
hyperAcc_bed <- hyperAcc[,c("Chr","Start","End","nucID","nucMACC","Strand")]
write.table(hyperAcc_bed , file="hyperAcc_monoNucs.bed",
            row.names = FALSE, sep="\t", quote=FALSE, col.names = FALSE)


#### full list containing all called nucleosomes
nucStats.full<-cbind(readCounts[,c("Chr", "Start","End", "Strand","GC_cont", "nucID")],
                     df[match(readCounts$nucID,df$nucID),c("nucMACC","category")],
                     regr.results[match(readCounts$nucID,rownames(regr.results)),])

write.table(nucStats.full, file="nucMACC_result_table_allNucs.tsv",
            row.names = FALSE, sep="\t", quote=FALSE,col.names = T)
