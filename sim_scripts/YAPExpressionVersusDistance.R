source("scripts/Definitions.R");

embryoData <- allData[[embryoHeading]];
embryos <- unique(embryoData);

treatments <- c(CONTROL, TREATED);

averageYAPRatioVersusDistance <- buildDataFrame(matrix(0.0,nrow=11,ncol=2),
                                                        11,
                                                        2,
                                                        CELL,
                                                        NULL,
                                                        c(paste(treatments[1], YAP_VERSUS_DISTANCE, sep="_"),
                                                          paste(treatments[2], YAP_VERSUS_DISTANCE, sep="_")));

colHeadings <- colnames(averageYAPRatioVersusDistance);
targetCol1 <- paste(NUCLEUS, YAP, MEAN_INTENSITY, sep ="_" );
targetCol2 <- paste(CYTOPLASM, YAP, MEAN_INTENSITY, sep ="_" );

for(e in embryos){
  thisData <- subset(allData, embryoData==e);
  n <- nrow(thisData);
  t <- thisData[[treatmentHeading]][1] + 1;
  currentYAPRatioVersusDistance <- rep(0.0, 11);
  
  counts <- rep(0, 11);
  for(i in 1:n){
    d <- round(thisData[[DISTANCE_TO_CENTRE]][i] * 10) + 1;
    counts[d] <- counts[d] + 1;
    currentYAPRatioVersusDistance[d] <- currentYAPRatioVersusDistance[d] + (thisData[[targetCol1]][i] / thisData[[targetCol2]][i]);
  }
  for(i in 1:length(counts)){
    if(counts[i] > 0){
      averageYAPRatioVersusDistance[[colHeadings[t]]][i] <- averageYAPRatioVersusDistance[[colHeadings[t]]][i] + currentYAPRatioVersusDistance[i] / counts[i];
    }
  }
}

write.csv(averageEmbryoMorphVersusDistance, paste("outputs", "gata3_versus_distance.csv", sep="/"));