source("scripts/Definitions.R");

outputHeadings <- c(TREATED, CONTROL);

embryoData <- allData[[embryoHeading]];
embryos <- unique(embryoData);
cellCounts <- list(data.frame(x=NULL), data.frame(x=NULL));

for(e in embryos){
  thisData <- subset(allData, embryoData==e);
  treatmentData <- thisData[[treatmentHeading]];
  row <- data.frame(matrix(ncol = 1, nrow = 0));
  colnames(row) <- outputHeadings[treatmentData[1] + 1];
  row[1,1] <- nrow(thisData);
  cellCounts[[treatmentData[1] + 1]] <- rbind(cellCounts[[treatmentData[1] + 1]], row);
}

cellCounts <- merge(cellCounts[[1]], cellCounts[[2]], by = 0, all = TRUE)[-1];

write.csv(cellCounts, paste("outputs", "cell_counts.csv", sep="/"));