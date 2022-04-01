#!/usr/bin/env Rscript

#### PACKAGE LOADING ####
library(RColorBrewer)

######## DATA LOADING ###########
args   <- commandArgs(TRUE)


heat.val<-read.delim(file = args[1])

#defining color palette
dark2 <- c(RColorBrewer::brewer.pal(8, "Dark2"),RColorBrewer::brewer.pal(8, "Set1"),RColorBrewer::brewer.pal(8, "Set2"))


#split label names 
lab<-strsplit(as.character(heat.val$bin.labels), split = "_")

#get conditions from profile

conditions<-NULL
names<-NULL
for(i in 2:length(lab)){
  if(lab[[i]][1]=="pooled"){
    next
  }
temp<-(lab[[i]][which(grepl( "U", lab[[i]], fixed = TRUE))])
assign(paste("data.",temp, sep=""), grep(paste(temp),heat.val$bin.labels, value = T))
names<-c(names,grep(paste(temp),heat.val$bin.labels, value = T))
conditions<-c(conditions,paste("data.",temp, sep=""))
}

names <- gsub("_monoNucs_profile","",names)


row.names(heat.val)<-as.character(heat.val$bin.labels)

TSS.pos<-which(colnames(heat.val)=="tick")

start.plot<-c(-1250)
end.plot<-c(1500)
pos<-seq(start.plot,end.plot,10)

profile.index<- (TSS.pos+start.plot/10):(TSS.pos+end.plot/10)

#get min value in 
for(i in 1:length(conditions)){
  curr_min<-min(heat.val[get(conditions[i]),profile.index ])
  if(i==1){
    low<-curr_min
  }
  else{
    if(curr_min<low){
      low<-curr_min
    }
  }
}

#get max value in 
for(i in 1:length(conditions)){
  curr_max<-max(heat.val[get(conditions[i]),profile.index ])
  if(i==1){
    high<-curr_max
  }
  else{
    if(curr_max>high){
      high<-curr_max
    }
  }
}


pdf("profile_monoNucs.pdf", width = 7, height=4)
par(mar=c(5.1, 4.1, 4.1, 11.8))
for(i in 1:length(conditions)){
  if(i==1){
    plot(pos,
      heat.val[get(conditions[i]),profile.index ], 
      type="l", 
      xlab="distance from TSS", ylab = "MNase fragment density",
      lwd=2, main="mono-nucs", col=dark2[i], ylim=c(low,high))
  }
  else{
    lines(pos,
      heat.val[get(conditions[i]), profile.index],
      col=dark2[i],lwd=2)
  }
}
  legend("topright",bty="n",
         legend=names, 
         col=dark2[1:length(conditions)], lwd=2,
       inset=c(-0.63,0), xpd=TRUE)
dev.off()


