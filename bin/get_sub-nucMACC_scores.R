#!/usr/bin/env Rscript

#### PACKAGE LOADING ####
library(LSD)
library(rtracklayer)

############################################
######## DATA LOADING ###########
args   <- commandArgs(TRUE)

#read sample input file
input<-read.csv(file=args[1])

## read count table
readCounts <- read.delim(file = args[2],header = T)
readCounts.Nucs<-read.delim(file = args[5],header = T)

## get lowest condition
lowest.cond<-args[3]
nucStats<-read.delim(args[4])

#get total counts statistics
featureCounts.mono<-read.delim(file=args[6])
colnames(featureCounts.mono)<-gsub("_mono.bam", "", colnames(featureCounts.mono))

featureCounts.sub<-read.delim(file = args[7])
colnames(featureCounts.sub)<-gsub("_sub.bam", "", colnames(featureCounts.sub))


dir.create("Figures", showWarnings = FALSE)


##############################################

#################################################
################### sub-nucMACC #################
#################################################
### DATA PREPARATION

#get count table
count.table<-readCounts[,c(6:(ncol(readCounts)-2))]
rownames(count.table)<-readCounts$nucID
colnames(count.table)<-gsub("_sub.bam", "", colnames(count.table))

### READ filtering based on raw reads in lowest concentration
#raw reads threshold
raw.flt<-5

#filter
idx.raw<-count.table[,lowest.cond]>=raw.flt
real.subNucs <- count.table[idx.raw,]

### READ COUNT NORMALIZATION

#get normalization factors for CPMs
normFactor<-apply(featureCounts.sub[,c(-1)],2,sum)/1e6

#Counts per million (normalization based on library size)
normCounts<-t(t(real.subNucs)/normFactor[match(colnames(real.subNucs),
                                               names(normFactor))]) 

### REGRESSION

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
    slope <- fit$coefficients[2]
    R2<-summary(fit)$r.squared
    return(data.frame(slope,R2))
}

regr.list <- apply(normCounts, 1, makeRegr)    #Actual calculation of linear regression
regr.results <-do.call(rbind, regr.list)

#####################################################
############# sub-nucMACC GC correction  ############
#####################################################

#extract sub-nucMACC
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
    heatscatter(x=gc.filt[ord],y=slope.filt[ord], cor=T, cexplot=0.5,
                ylab="sub-nucMACC scores",
                xlab="GC content")
    abline(h=0,col="black",lwd=3,lty=2)
    lines(gc.filt[ord],predict.slope[ord], col="#5CB85B",lwd=3)
dev.off()

#### PLOT OF CORRELATION AFTER THE GC NORMALIZATION
png("Figures/after_gc_norm_finalPlot.png", width = 1280, height = 1280, res =300 )
    heatscatter(x=gc.filt[ord], y=(slope.filt[ord])-(predict.slope[ord]-median(slope.filt)),
                cor=T,cexplot=0.5, ylab="sub-nucMACC scores",
                xlab="GC content")
    abline(h=0,col="black",lwd=3,lty=2)
dev.off()


# GC normalize sub-nucMACC scores
sub_nucMACC_scores <- slope.filt-(predict.slope-median(slope.filt))


########################################
#### create output list and bedgraph ###
########################################

subnucStats.full<-cbind(readCounts[,c("Chr", "Start","End", "Strand","GC_cont", "nucID")],
                     "sub.nucMACC"=sub_nucMACC_scores[match(readCounts$nucID,
                                                            names(sub_nucMACC_scores))],
                     regr.results[match(readCounts$nucID,rownames(regr.results)),])

#list filtered by raw read count and GC
subnucStats<-subnucStats.full[!is.na(subnucStats.full$sub.nucMACC),]

#export sub-nucMACC as bedgraph
bedgraph <- subnucStats[,c("Chr","Start","End","sub.nucMACC")]
write.table(bedgraph , file="sub-nucMACC_scores.bedgraph",
            row.names = FALSE, sep="\t", quote=FALSE, col.names = FALSE)




#########################################
###### comparison with nucMACC scores ###
#########################################

#convert to GRanges
subnuc.gr<-GRanges(subnucStats)
nuc.gr<-GRanges(nucStats)

# get unique subnuc positions
subNucs_unique  <- subsetByOverlaps(subnuc.gr,nuc.gr,
                                    minoverlap = 70,invert=TRUE)


###################################################################
############### Finding positions enriched in subnucs##############
###################################################################

## overlapping
ovrlp<-findOverlaps(subnuc.gr,nuc.gr, minoverlap = 70)
sub.nucMACC.ovrl<-subnuc.gr[ovrlp@from]
nucMACC.ovrl<-nuc.gr[ovrlp@to]

#get count table of monoNucs
count.table.Nucs<-readCounts.Nucs[,c(6:(ncol(readCounts.Nucs)-2))]
rownames(count.table.Nucs)<-readCounts.Nucs$nucID
colnames(count.table.Nucs)<-gsub("_mono.bam", "", colnames(count.table.Nucs))

##############################################################
########### norm by total filtered fragments per million #####

sub.total<-sum(featureCounts.sub[,lowest.cond])/1e6
nuc.total<-sum(featureCounts.mono[,lowest.cond])/1e6

cpms.nuc<-count.table.Nucs[as.character(nucMACC.ovrl$nucID),lowest.cond]/nuc.total
names(cpms.nuc)<-as.character(nucMACC.ovrl$nucID)
cpms.subnuc<-count.table[as.character(sub.nucMACC.ovrl$nucID),lowest.cond]/sub.total
names(cpms.subnuc)<-as.character(sub.nucMACC.ovrl$nucID)

#get enrichment in
ratio<-log2(cpms.subnuc+0.01)-log2(cpms.nuc+0.01)

## select more than 4 fold higher coverage in lowest MNase condition (log2)
enriched.sub<-names(ratio)[ratio>2]


########################################################
########## get fragile- and non-canonical Nucs #########
########################################################
#filter relevant subnucleosomes
selected.subNucs<-c(enriched.sub,as.character(subNucs_unique$nucID)) # subnucs selected
sub_nucMACC_scores.select<-sub_nucMACC_scores #[selected.subNucs]

# prepare for calling
# normalize nucMACC scores to maximum 1
sub_nucMACC_sort <- sort(sub_nucMACC_scores.select)
sub_nucMACC_norm <- (sub_nucMACC_sort/(max(sub_nucMACC_sort)-min(sub_nucMACC_sort)))
#normalize ranks to maximum 1
sub_nucMACC_rank <- rank(sub_nucMACC_sort, ties.method = "first")/length(sub_nucMACC_sort)

### get LOESS smoothing function
#define to take 1000 nucleosomes as window for loess
span.loess.pre<-1000/length(sub_nucMACC_scores)
span.loess<-ifelse(span.loess.pre<0.05, span.loess.pre, 0.05)
lm <- loess(sub_nucMACC_norm~sub_nucMACC_rank,span=span.loess)

#Curve fitting
loess.line <- predict(lm)

#get first derivative to deduce the slope of loess function
loess.f1 <- diff(loess.line)/diff(sub_nucMACC_rank)

#set the slope cutOff to 1
cutOff1 <- min(which(loess.f1 < 1))
cutOff2 <- max(which(loess.f1 < 1))

## show loess derivative and cutoff 1
png("Figures/loess_first_order_derivative.png",
        width = 1280, height = 1280, res =300 )
    plot(loess.f1, xlab="sub-nucMACC rank", ylab="loess first order derivative",
         type="l")
    abline(h = 1, col="#D43F39")
dev.off()

#final selection
png("Figures/nucMACC_selection.png",
    width = 1280, height = 1280, res =300 )
    plot(sub_nucMACC_rank,sub_nucMACC_norm, ylab="normalized sub-nucMACC score",
         xlab = "normalized rank", bty="n", pch=15)
    abline(v=sub_nucMACC_rank[cutOff1], col="#2166AC", lty=2)
    abline(v=sub_nucMACC_rank[cutOff2], col="#B2182B", lty=2)
    points(sub_nucMACC_rank[1:cutOff1], sub_nucMACC_norm[1:cutOff1],pch=15,
           col="#2166AC")
    points(sub_nucMACC_rank[cutOff2:length(sub_nucMACC_rank)],
           sub_nucMACC_norm[cutOff2:length(sub_nucMACC_rank)],pch=15,col="#B2182B")
    lines(sub_nucMACC_rank,loess.line, col="#D95F02")
dev.off()



### select the nucs with nucMACC scores deviating from the mean
#scores need to be positive or negative 
cutoff.low<-min(c(sub_nucMACC_sort[cutOff1],0))
cutoff.high<-max(c(sub_nucMACC_sort[cutOff2],0))

sub_nucMACC_low.pre <- sub_nucMACC_scores.select < cutoff.low
sub_nucMACC_high.pre <- sub_nucMACC_scores.select > cutoff.high


sub_nucMACC_low <- sub_nucMACC_low.pre & (names(sub_nucMACC_low.pre) %in% selected.subNucs)
sub_nucMACC_high <- sub_nucMACC_high.pre & (names(sub_nucMACC_high.pre) %in% selected.subNucs)

### generate data table
df<-data.frame(nucID=names(sub_nucMACC_scores), sub.nucMACC=sub_nucMACC_scores)

# define special
df$category<-"not selected"
df[names(sub_nucMACC_scores.select),]$category<-"normal"

df[names(sub_nucMACC_scores.select),]$category[sub_nucMACC_low]<-"un-stable"
df[names(sub_nucMACC_scores.select),]$category[sub_nucMACC_high]<-"non-canonical"

#get selection criteria
df$selection<-"not selected"
df[enriched.sub,"selection"] <- "enriched"
df[as.character(subNucs_unique$nucID),"selection"] <- "unique"

##concatenate tables
subnucStats<-data.frame(subnucStats,df[,c("category","selection")])

write.table(subnucStats, file="sub-nucMACC_result_table.tsv",
            row.names = FALSE, sep="\t", quote=FALSE,col.names = T)


#export bed file and define score nucMACC
#non-canonical
noCan<-subset(subnucStats, category=="non-canonical")
noCan_bed <- noCan[,c("Chr","Start","End","nucID","sub.nucMACC","Strand")]
write.table(noCan_bed , file="nonCanonical_subNucs.bed",
            row.names = FALSE, sep="\t", quote=FALSE, col.names = FALSE)

#un-stable
unStab<-subset(subnucStats, category=="un-stable")
unStab_bed <- unStab[,c("Chr","Start","End","nucID","sub.nucMACC","Strand")]
write.table(unStab_bed , file="unStable_subNucs.bed",
            row.names = FALSE, sep="\t", quote=FALSE, col.names = FALSE)
