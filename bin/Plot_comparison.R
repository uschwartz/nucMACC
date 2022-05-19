#!/usr/bin/env Rscript

### Load Packages 
library(ggplot2)

# Reading the output-files from GenerateTxtFragCounts.R 
txt_output <- list.files(path = ".", pattern = "fragment_statistic.txt", full.names = TRUE)

# Generate a matrix that contains all the information to be plotted. 
# Here the absolute counts are needed to calculate the steps lost in each step of the pipeline.
lost_matrix <- matrix(data = NA, nrow = 6, ncol = 2)
names <- c()
for (i in 1:length(txt_output)){
  df <- read.table(txt_output[i], header = TRUE, sep = " ")
  lost_matrix[1, i] <- df[1,1]
  lost_matrix[2, i] <- df[1, 1] - df[1,2]
  lost_matrix[3, i] <- (df[1, 1] - df[1,3]) - lost_matrix[2, i]
  lost_matrix[4, i] <- (df[1, 1] - (df[1, 4] + df[1,5])) - (lost_matrix[3, i] + lost_matrix[2, i])
  lost_matrix[5, i] <- df[1, 4]
  lost_matrix[6, i] <- df[1, 5]
  rownames(lost_matrix) <- c("Sequenced", "Not aligned", "Quality-filtered", "Size- and Blacklist-filtered", colnames(df[, 4:5]))
  names[i] <- rownames(df)
}
colnames(lost_matrix) <- names
destination <- c()
out <- c()
for (i in 1:ncol(lost_matrix)){
  destination <- c(destination, rep(colnames(lost_matrix)[i], nrow(lost_matrix)))
}

# Further editing of the matrix to a data frame with all specifications needed to plot using ggplot2.
fragments <- c(lost_matrix[, 1], lost_matrix[,2])
plot_matrix <- cbind(fragments, destination)
plot_matrix <- cbind(c(rownames(lost_matrix), rownames(lost_matrix)), plot_matrix)
colnames(plot_matrix) <- c("Type", "Counts", "Sample")
rownames(plot_matrix) <- c(1:nrow(plot_matrix))
plot_frame <- as.data.frame(plot_matrix)
to_keep <- ifelse(plot_frame[,1] != "Sequenced", FALSE, TRUE)
plot_frame <- plot_frame[!to_keep, ]
plot_frame$Type <- factor(plot_frame$Type, levels = c("Not aligned", "Quality-filtered",  "Size- and Blacklist-filtered", "MonoNuc",  "SubNuc"))
plot_frame$Counts <- as.numeric(plot_frame$Counts)

g <- ggplot(plot_frame, aes(x = Sample, y = Counts/10**6, fill = Type)) + 
      geom_bar(stat = "identity") +
      scale_fill_manual(values=c("gray71", "gray60", "gray40", "coral", "coral3")) + 
      scale_y_continuous(breaks = seq(0, max(lost_matrix)/10**6, by = 50)) +
      xlab("") +
      ylab("read fragments per million") +
      labs(fill = "") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.y = element_text(color = "grey30", size = 10))
pdf("read_statistic.pdf", width = ncol(lost_matrix) + 2, height = 4)
print(g)
dev.off()
