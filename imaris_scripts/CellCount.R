source("sim_scripts/Definitions.R");

pdf(paste("plots", "centroid_error_versus_snr.pdf", sep=.Platform$file.sep));
boxplot(Centroid_Error~snr, data=allData, main="Centroid Error", xlab="SNR", ylab="Centroid Error (Microns)", notch=TRUE, col=(c("red", "green", "blue")));

dev.off();

cellCountErrors <- data.frame();

snrs <- unique(allData[[SNR]]);

for(i in snrs){
  temp<-allData[allData$snr == i, ];
  temp<-temp[!duplicated(temp$Index),];
  cellCountErrors <- rbind(cellCountErrors, data.frame(SNR=i,
                                                       Mean_Cell_Count_Error=mean(temp[, CELL_COUNT_ERROR]),
                                                       SD_Cell_Count_Error=sd(temp[, CELL_COUNT_ERROR])));
}

write.csv(cellCountErrors, paste("outputs", "cell_count_errors.csv", sep=.Platform$file.sep));