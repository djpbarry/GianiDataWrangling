source("scripts/Definitions.R");

outputHeadings <- c(CONTROL, TREATED);

embryoData <- allData[[embryoHeading]];
embryos <- unique(embryoData);
cellCounts <- list(data.frame(x=NULL), data.frame(x=NULL));

for(e in embryos){
  thisData <- subset(allData, embryoData==e);
  treatmentData <- thisData[[treatmentHeading]];
  row <- data.frame(matrix(ncol = 1, nrow = 0));
  colnames(row) <- outputHeadings[treatmentData[1]];
  row[1,1] <- nrow(thisData);
  cellCounts[[treatmentData[1]]] <- rbind(cellCounts[[treatmentData[1]]], row);
}

cellCounts <- merge(cellCounts[[CONTROL_VALUE]], cellCounts[[TREATED_VALUE]], by = 0, all = TRUE)[-1];

write.csv(cellCounts, paste("outputs", "cell_counts.csv", sep="/"));

saveConfidenceIntervals(cellCounts[[CONTROL_VALUE]][!is.na(cellCounts[[CONTROL_VALUE]])], cellCounts[[TREATED_VALUE]][!is.na(cellCounts[[TREATED_VALUE]])], "outputs", "cell_counts_CI.csv");

pdf(paste("plots", "cell_counts.pdf", sep=.Platform$file.sep));

hist(cellCounts[[CONTROL_VALUE]][!is.na(cellCounts[[CONTROL_VALUE]])], col=rgb(0.8,0.0,0.0,0.5), freq = TRUE, border="black", main="Cell Counts", ylab="Number of Embryos",xlab="Number of Cells", xlim=c(10,40), ylim=c(0,10), breaks=seq(from = 10, to = 40, by = 2.5));
hist(cellCounts[[TREATED_VALUE]][!is.na(cellCounts[[TREATED_VALUE]])], col=rgb(0.0,0.8,0.0,0.5), freq = TRUE, border="black", breaks=seq(from = 10, to = 40, by = 2.5), add=TRUE);
legend("topright", c(CONTROL, TREATED), fill=c("red", "green"));

dev.off();