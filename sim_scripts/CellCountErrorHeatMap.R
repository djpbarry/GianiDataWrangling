source("sim_scripts/Definitions.R");

snrs <- unique(allData[[SNR]]);
cellCounts <- unique(allData[[GROUND_TRUTH_N]]);

snrs<-sort(snrs);
cellCounts<-sort(cellCounts);

cellCountErrors <- data.frame();

for(i in snrs){
  temp1<-allData[allData[[SNR]] == i, ];
  for(j in cellCounts){
    temp2 <- temp1[temp1$Ground_Truth_Cell_Count == j, ];
    errors <- temp2[!is.nan(temp2$Centroid_Error),];
    cellCountErrors <- rbind(cellCountErrors, c(i, j, mean(errors$Centroid_Error), mean(abs(errors$Normalised_Volume_Error)), mean(abs(errors$Proportionate_Volume_Error))));
  }
}

colnames(cellCountErrors) <- c(SNR, "N_CELLS", "MEAN_CENTROID_ERROR", "VOLUME_ERROR", "CORRECTED_VOLUME_ERROR");

#pdf(paste("plots", "sim_cell_count_errors.pdf", sep=.Platform$file.sep));

#heatmap(as.matrix(cellCountErrors), Colv = NA, Rowv = NA, col=colorRampPalette(brewer.pal(8, "RdYlBu"))(16), xlab="SNR", ylab="Cell Number", cexCol=1.0);

#dev.off();

heatMapData <- allData[!duplicated(allData[[INDEX]]),];
axislabel <- element_text(hjust=0.5, size=18);

p <- ggplot(heatMapData, aes(as.character(snr), Ground_Truth_Cell_Count, fill=Proportional_Cell_Count_Error)) + 
  geom_tile() + scale_fill_gradient(low="blue", high="red") + xlab("SNR") + ylab("Number of Cells") + labs(fill = "Cell Count Error");
p + theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel);

p <- ggplot(cellCountErrors, aes(as.character(snr), N_CELLS, fill=MEAN_CENTROID_ERROR)) + 
  geom_tile() + scale_fill_gradient(low="blue", high="red") + xlab("SNR") + ylab("Number of Cells") + labs(fill = "Localisation Error");
p + theme(axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel);

p <- ggplot(cellCountErrors, aes(as.character(snr), N_CELLS, fill=VOLUME_ERROR)) + 
  geom_tile() + scale_fill_gradient(low="blue", high="red") + xlab("SNR") + ylab("Number of Cells") + labs(fill = "Segmentation Error");
p + theme(axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel);

p <- ggplot(cellCountErrors, aes(as.character(snr), N_CELLS, fill=CORRECTED_VOLUME_ERROR)) + 
  geom_tile() + scale_fill_gradient(low="blue", high="red") + xlab("SNR") + ylab("Number of Cells") + labs(fill = "Segmentation Error");
p + theme(axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel);

#plot(x, controlData, xlab=xLabel, main=title, ylab=yLabel, ylim=yLimits, col="red", pch=15);
plot(allData$Normalised_Distance_to_Centre, allData$Normalised_Volume_Error);
