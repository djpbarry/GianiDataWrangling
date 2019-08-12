source("scripts/Definitions.R");

embryoData <- allData[[embryoHeading]];
embryos <- unique(embryoData);

treatments <- c(CONTROL, TREATED);

averageGATAVersusDistance <- buildDataFrame(matrix(0.0,nrow=11,ncol=2),
                                                        11,
                                                        2,
                                                        NUCLEUS,
                                                        NULL,
                                                        c(paste(treatments[1], GATA_VERSUS_DISTANCE, sep="_"),
                                                          paste(treatments[2], GATA_VERSUS_DISTANCE, sep="_")));

colHeadings <- colnames(averageGATAVersusDistance);
targetCol <- paste(NUCLEUS, GATA, MEAN_INTENSITY, sep ="_" );

for(e in embryos){
  thisData <- subset(allData, embryoData==e);
  n <- nrow(thisData);
  t <- thisData[[treatmentHeading]][1] + 1;
  currentGATAVersusDistance <- rep(0.0, 11);
  
  counts <- rep(0, 11);
  for(i in 1:n){
    d <- round(thisData[[DISTANCE_TO_CENTRE]][i] * 10) + 1;
    counts[d] <- counts[d] + 1;
    currentGATAVersusDistance[d] <- currentGATAVersusDistance[d] + thisData[[targetCol]][i];
  }
  for(i in 1:length(counts)){
    if(counts[i] > 0){
      averageGATAVersusDistance[[colHeadings[t]]][i] <- averageGATAVersusDistance[[colHeadings[t]]][i] + currentGATAVersusDistance[i] / counts[i];
    }
  }
}

write.csv(averageEmbryoMorphVersusDistance, paste("outputs", "gata3_versus_distance.csv", sep="/"));