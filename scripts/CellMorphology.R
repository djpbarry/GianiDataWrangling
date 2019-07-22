source("scripts/Definitions.R");

embryoData <- allData[[embryoHeading]];
embryos <- unique(embryoData);
cellMorph <- list(data.frame(x=NULL), data.frame(x=NULL));

for(e in embryos){
  thisData <- subset(allData, embryoData==e);
  treatmentData <- thisData[[treatmentHeading]];
  if(treatmentData[1] > 0){
    label = CONTROL;
  } else {
    label = TREATED;
  }
  row <- morphData <- buildDataFrame(c(mean(thisData[[NUCLEUS_SURFACE_AREA_MIC]]),
                                       mean(thisData[[NUCLEUS_VOLUME_MIC]]),
                                       mean(thisData[[CELL_SURFACE_AREA_MIC]]),
                                       mean(thisData[[CELL_VOLUME_MIC]])),
                                     1,
                                     4,
                                     label,
                                     NULL,
                                     c(paste("Mean", NUCLEUS_SURFACE_AREA_MIC, sep="_"),
                                       paste("Mean", NUCLEUS_VOLUME_MIC, sep="_"),
                                       paste("Mean", CELL_SURFACE_AREA_MIC, sep="_"),
                                       paste("Mean", CELL_VOLUME_MIC, sep="_")));
  cellMorph[[treatmentData[1] + 1]] <- rbind(cellMorph[[treatmentData[1] + 1]], row);
}

cellMorph <- merge(cellMorph[[1]], cellMorph[[2]], by = 0, all = TRUE)[-1];

write.csv(cellMorph, paste("outputs", "cell_morphology.csv", sep="/"));