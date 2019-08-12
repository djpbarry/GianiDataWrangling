source("scripts/Definitions.R");

embryoData <- allData[[embryoHeading]];
embryos <- unique(embryoData);
cellMorph <- list(data.frame(x=NULL), data.frame(x=NULL));
NUM_EMB <- "Number of Embryos";

for(e in embryos){
  thisData <- subset(allData, embryoData==e);
  treatmentData <- thisData[[treatmentHeading]];
  if(treatmentData[1] > CONTROL_VALUE){
    label = TREATED;
  } else {
    label = CONTROL;
  }
  row <- morphData <- buildDataFrame(c(mean(thisData[[NUCLEUS_VOLUME_MIC]]),
                                       mean(thisData[[CELL_VOLUME_MIC]]),
                                       mean(thisData[[NUCLEUS_VOLUME_MIC]] / thisData[[CELL_VOLUME_MIC]])),
                                     1,
                                     3,
                                     label,
                                     NULL, c(paste(MEAN, NUCLEUS_VOLUME_MIC, sep="_"),
                                             paste(MEAN, CELL_VOLUME_MIC, sep="_"),
                                             paste(MEAN, NUC_TO_CELL_VOLUME_RATIO, sep="_")));
  cellMorph[[treatmentData[1]]] <- rbind(cellMorph[[treatmentData[1]]], row);
}

columns <- c(colnames(cellMorph[[CONTROL_VALUE]]), colnames(cellMorph[[TREATED_VALUE]]));

lightRed <- rgb(0.8,0.0,0.0,0.5);
lightGreen <- rgb(0.0,0.8,0.0,0.5);

saveHistogram("plots", "nuc_vol.pdf", cellMorph[[CONTROL_VALUE]][[columns[1]]], cellMorph[[TREATED_VALUE]][[columns[4]]], lightRed, lightGreen, "Nuclear Volume", expression(paste("Mean Nuclear Volume (",mu,"m^3)")), NUM_EMB, c(1000,3500), c(0,10), seq(from=1000, to=3500, by=250));
saveHistogram("plots", "cell_vol.pdf", cellMorph[[CONTROL_VALUE]][[columns[2]]], cellMorph[[TREATED_VALUE]][[columns[5]]], lightRed, lightGreen, "Cell Volume", expression(paste("Mean Cell Volume (",mu,"m^3)")), NUM_EMB, c(4000,18000), c(0,8), seq(from=4000, to=18000, by=2000));
saveHistogram("plots", "nuc_to_cell_vol_ratio.pdf", cellMorph[[CONTROL_VALUE]][[columns[3]]], cellMorph[[TREATED_VALUE]][[columns[6]]], lightRed, lightGreen, "Nuclear/Cell Volume Ratio", "Mean Nuclear:Cell Volume Ratio", NUM_EMB, c(0.1,0.5), c(0,10), seq(from=0.1, to=0.5, by=0.05));

saveConfidenceIntervals(cellMorph[[CONTROL_VALUE]][[columns[1]]][!is.na(cellMorph[[CONTROL_VALUE]][[columns[1]]])], cellMorph[[TREATED_VALUE]][[columns[4]]][!is.na(cellMorph[[TREATED_VALUE]][[columns[4]]])], "outputs", "nuc_vol_CI.csv");
saveConfidenceIntervals(cellMorph[[CONTROL_VALUE]][[columns[2]]][!is.na(cellMorph[[CONTROL_VALUE]][[columns[2]]])], cellMorph[[TREATED_VALUE]][[columns[5]]][!is.na(cellMorph[[TREATED_VALUE]][[columns[5]]])], "outputs", "cell_vol_CI.csv");
saveConfidenceIntervals(cellMorph[[CONTROL_VALUE]][[columns[3]]][!is.na(cellMorph[[CONTROL_VALUE]][[columns[3]]])], cellMorph[[TREATED_VALUE]][[columns[6]]][!is.na(cellMorph[[TREATED_VALUE]][[columns[6]]])], "outputs", "nuc_to_cell_vol_CI.csv");

cellMorph <- merge(cellMorph[[CONTROL_VALUE]], cellMorph[[TREATED_VALUE]], by = 0, all = TRUE)[-1];

write.csv(cellMorph, paste("outputs", "cell_morphology.csv", sep="/"));
