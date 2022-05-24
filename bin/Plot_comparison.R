#!/usr/bin/env Rscript

### Load Packages 
library(ggplot2)
library(reshape2)

# Reading the output-files from GenerateTxtFragCounts.R 
txt_output <- list.files(path = ".", pattern = "fragment_statistic.txt", full.names = T)

# Generate a matrix that contains all the information to be plotted. 
# Here the absolute counts are needed to calculate the steps lost in each step of the pipeline.
lost_matrix <- matrix(data = NA, nrow = 6, ncol = length(txt_output))
rownames(lost_matrix) <- c("Sequenced", "Not aligned", "Quality-filtered",
                           "Size- and Blacklist-filtered", "SubNuc","MonoNuc")

names <- c()
txt_summary <- matrix(data = NA, nrow = length(txt_output), ncol = 5)
for (i in 1:length(txt_output)){
  df <- read.table(txt_output[i], header = T, sep = "\t")
  lost_matrix[1, i] <- df$Sequenced
  lost_matrix[2, i] <- df$Sequenced - df$Aligned
  lost_matrix[3, i] <- (df$Sequenced - df$MAPQC.filtered) -
    lost_matrix["Not aligned", i]
  lost_matrix[4, i] <- (df$Sequenced - (df$SubNuc + df$MonoNuc)) - 
    (lost_matrix["Quality-filtered", i] + lost_matrix["Not aligned", i])
  lost_matrix[5, i] <- df$SubNuc
  lost_matrix[6, i] <- df$MonoNuc
  names[i] <- rownames(df)
  colnames(txt_summary) <- colnames(df)
  txt_summary[i, ] <- unlist(df)
}
colnames(lost_matrix) <- names
rownames(txt_summary) <- names

### Plot the relative amount of reads based on the originally sequenced reads. 
relative <- matrix(data = NA, nrow = nrow(lost_matrix), ncol = ncol(lost_matrix))
for (i in 1:ncol(lost_matrix)){
  relative[,i] <- sapply(lost_matrix[,i], function(x){(x/lost_matrix[1,i])*100})
}
rownames(relative) <- rownames(lost_matrix)
colnames(relative) <- colnames(lost_matrix)
rel_frame <- melt(relative[-1,], varnames = c("Type", "Sample"), value.name = "Counts")
rel_frame$Type <- factor(rel_frame$Type, levels = c("Not aligned", "Quality-filtered",  
                                                    "Size- and Blacklist-filtered", "MonoNuc",  "SubNuc"))
rel_frame$Counts <- as.numeric(rel_frame$Counts)
r <- ggplot(rel_frame, aes(x = Sample, y = Counts, fill = Type)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("gray71", "gray60", "gray40", "coral", "coral3")) + 
  xlab("") +
  ylab("% per fragments sequenced") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.y = element_text(color = "grey30", size = 10))

### Plot the absolute amounts of reads. 
plot_frame <- melt(lost_matrix[-1,], varnames = c("Type", "Sample"), value.name = "Counts")
plot_frame$Type <- factor(plot_frame$Type, levels = c("Not aligned", "Quality-filtered",  
                                                      "Size- and Blacklist-filtered", "MonoNuc",  "SubNuc"))
plot_frame$Counts <- as.numeric(plot_frame$Counts)
g <- ggplot(plot_frame, aes(x = Sample, y = Counts/10**6, fill = Type)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("gray71", "gray60", "gray40", "coral", "coral3")) + 
  scale_y_continuous(breaks = seq(0, max(lost_matrix)/10**6, by = 50)) +
  xlab("") +
  ylab("fragments per million") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.y = element_text(color = "grey30", size = 10))

### Save both plots to a summary pdf.
pdf("fragment_statistic.pdf", width = ncol(lost_matrix) + 2, height = 4)
print(g)
print(r)
dev.off()

### Save a further summary file in tab-delimited txt-format to give statistics overview.
txt_summary.ext<-cbind("Samples"=rownames(txt_summary), txt_summary)
write.table(x = txt_summary.ext, file = "Fragment_statistic_summary.txt",
            append = F, sep = "\t", dec = ".", row.names = F, 
            col.names = T, quote = F)
