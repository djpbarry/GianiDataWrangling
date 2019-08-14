source("sim_scripts/Definitions.R");

#write.csv(cellCounts, paste("outputs", "cell_counts.csv", sep="/"));

#saveConfidenceInterval(cbind(cellCounts[[CELL_COUNT_ERROR]][cellCounts[[SNR]] == 1], cellCounts[[CELL_COUNT_ERROR]][cellCounts[[SNR]] == 2], cellCounts[[CELL_COUNT_ERROR]][cellCounts[[SNR]] == 3]), c("1", "2", "3"), "outputs", "cell_count_error_CI.csv");

pdf(paste("plots", "cell_count_error_versus_snr.pdf", sep=.Platform$file.sep));
boxplot(Cell_Count_Error~snr, data=allData[!duplicated(allData$Index),], main="Cell Count Error", xlab="SNR", ylab="Cell Count Error", notch=TRUE, col=(c("darkgreen")));
legend("topright", c("2"), fill=c("darkgreen"));

dev.off();

pdf(paste("plots", "centroid_error_versus_snr.pdf", sep=.Platform$file.sep));
boxplot(Centroid_Error~snr, data=allData, main="Centroid Error", xlab="SNR", ylab="Centroid Error (Microns)", notch=TRUE, col=(c("darkgreen")));
legend("topright", c("2"), fill=c("darkgreen"));

dev.off();

pdf(paste("plots", "cell_count_error_versus_cell_number.pdf", sep=.Platform$file.sep));
boxplot(Cell_Count_Error~Ground_Truth_Cell_Count, allData[!duplicated(allData$Index),], main="Cell Count Error", xlab="Cell Number", ylab="Cell Count Error", notch=TRUE, col=(c("darkgreen")));
legend("topright", c("2"), fill=c("darkgreen"));

dev.off();
