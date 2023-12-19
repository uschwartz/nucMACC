#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
###############################################################################  
#####################  get best positioned nucleosome reference ###############
###############################################################################  
library(rtracklayer)

input <- read.csv(args[1])
path_analysis <- unique(input$path_analysis)
bedfiles<-c(paste0(path_analysis[1], "/RUN/05_nucMACC/hyperAcc_monoNucs.bed"),
            paste0(path_analysis[2], "/RUN/05_nucMACC/hyperAcc_monoNucs.bed"),
            paste0(path_analysis[1], "/RUN/05_nucMACC/hypoAcc_monoNucs.bed"),
            paste0(path_analysis[2], "/RUN/05_nucMACC/hypoAcc_monoNucs.bed"))
## load bed files
bed.list<-GRangesList()
for (i in bedfiles){
bed.list[[which(bedfiles==i)]]<-import.bed(i)
}


### merge bed files together to one reference file; condition if two positions
### overlap > 100bp stitch them together

#start with first entry as query
newquery<-bed.list[[1]]
#remove from list
timepoint_peaks<-(bed.list[(-1)])


p_reduce<-GRanges()
# iterate over all timepoints
for(a in seq_along(timepoint_peaks)) { 
    print(c("This is round #", a))
  
    query1<-newquery
    subject1<-timepoint_peaks[[a]]
    
    
    # create Hits object, to later exclude them from query/ subject 
    temp_overlap <- findOverlaps(query1, subject1, minoverlap = 100)
    
    
    # create Pairs object to reduce regions pairwise
    temp_pairs <- findOverlapPairs(query1, subject1, minoverlap = 100)
    
    
    # get regions that are only in query/ subject to include them later
    if (isEmpty(queryHits(temp_overlap))==T){
      only_A<-query1
      only_B<-subject1
    }else{
      only_A<-query1[-queryHits(temp_overlap)]
      only_B<-subject1[-subjectHits(temp_overlap)]
    }
    
    ## test exact overlap of regions
    exactOvrl<-temp_pairs@first[(temp_pairs@first %in% temp_pairs@second)]
    temp_pairs2<-temp_pairs[!(temp_pairs@first %in% temp_pairs@second)]
    
    
    # reduce pairwise (get region from query.start to subject.end)
    start_temp <- start(temp_pairs2@first)
    start_second <- which(start(temp_pairs2@second) < start(temp_pairs2@first))
    
    start_temp[start_second]<- start(temp_pairs2@second)[start_second]
    
    end_temp <- end(ranges(temp_pairs2@first))
    end_second <- which(end(temp_pairs2@second) > end(temp_pairs2@first))
    
    end_temp[end_second] <- end(temp_pairs2@second)[end_second]
    
    s <- reduce(GRanges(seqnames = seqnames(temp_pairs2@first),
                        ranges = IRanges(start = (start_temp+50),
                                         end = (end_temp-50))))
    start(s)<-(start(s)-50)
    end(s)<-(end(s)+50)
    
    p_reduce <- c(p_reduce, s)
    
    # create query for next comparison step with regions: Overlapping, only in query, only in subject
    newquery <- c(p_reduce, exactOvrl ,only_A, only_B)
    # empty GRanges again for next pairwise reduction
    p_reduce <- GRanges()
        
}

## merge reference peaks if they share 120 bp
nucRef<-newquery
OvrlpCount<-countOverlaps(nucRef,nucRef, minoverlap = 120)
ovrlps<-which(OvrlpCount>1)

if(length(ovrlps)>0){
    merged<-reduce(nucRef[ovrlps])
    newquery<-c(nucRef[-ovrlps],merged)
}

start(newquery)<- start(newquery)-1
# write final result to bed file
base<-as.data.frame(newquery)
# therefore bring to bed file order
base_w<-data.frame(base[,1:3], base$name, base$score, base$strand)
write.table(base_w, file= "NucPosRef_allNucs.bed", 
            row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE)





