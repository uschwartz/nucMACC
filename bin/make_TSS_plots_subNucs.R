#!/usr/bin/env Rscript

#### PACKAGE LOADING ####
library(RColorBrewer)

######## DATA LOADING ###########
args   <- commandArgs(TRUE)

heat.val<-read.delim(file = args[1])
#defining color palette
dark2 <- c(RColorBrewer::brewer.pal(8, "Dark2"),RColorBrewer::brewer.pal(8, "Set1"),RColorBrewer::brewer.pal(8, "Set2"))


names <- gsub("_subNucs_profile","",as.character(heat.val$bin.labels)[-1])

row.names(heat.val)<-c("bin",names)

TSS.pos<-which(colnames(heat.val)=="tick")

start.plot<-c(-1250)
end.plot<-c(1500)
pos<-seq(start.plot,end.plot,10)

profile.index<- (TSS.pos+start.plot/10):(TSS.pos+end.plot/10)

#get min value in 
for(i in names){
  curr_min<-min(heat.val[i,profile.index ])
  curr_max<-max(heat.val[i,profile.index ])
  if(i==names[1]){
    low<-curr_min
    high<-curr_max
  }
  else{
      low<-min(curr_min,low)
      high<-max(curr_max,high)
  }
}



pdf("profile_subNucs.pdf", width = 7, height=4)
par(mar=c(5.1, 4.1, 4.1, 11.8))
for(i in 1:length(names)){
  if(i==1){
    plot(pos,
         heat.val[names[i],profile.index ], 
         type="l", 
         xlab="distance from TSS", ylab = "MNase fragment density",
         lwd=2, main="sub-nucs", col=dark2[i], ylim=c(low,high))
  }
  else{
    lines(pos,
          heat.val[names[i], profile.index],
          col=dark2[i],lwd=2)
  }
}
legend("topright",bty="n",
       legend=names, 
       col=dark2[1:length(names)], lwd=2,
       inset=c(-0.37,0), xpd=TRUE)
dev.off()


