source("scripts/Definitions.R");

embryoData <- allData[[embryoHeading]];
embryos <- unique(embryoData);

treatments <- c(CONTROL, TREATED);

averageEmbryoMorphVersusDistance <- buildDataFrame(matrix(0.0,nrow=11,ncol=2),
                                                        11,
                                                        2,
                                                        NULL,
                                                        NULL,
                                                        c(paste(treatments[1], MEAN, NUCLEAR_SURFACE_AREA_MICRONS_VERSUS_DISTANCE, sep="_"),
                                                          paste(treatments[2], MEAN, NUCLEAR_SURFACE_AREA_MICRONS_VERSUS_DISTANCE, sep="_")));

for(e in embryos){
  thisData <- subset(allData, embryoData==e);
  n <- nrow(thisData);
  t <- thisData[[treatmentHeading]][1] + 1;
  colHeading <- paste(treatments[t], MEAN, NUCLEAR_SURFACE_AREA_MICRONS_VERSUS_DISTANCE, sep="_");
  currentEmbryoMorphVersusDistance <- rep(0.0, 11);
  
  counts <- rep(0, 11);
  for(i in 1:n){
    d <- round(thisData[[DISTANCE_TO_CENTRE]][i] * 10) + 1;
    counts[d] <- counts[d] + 1;
    currentEmbryoMorphVersusDistance[d] <- currentEmbryoMorphVersusDistance[d] + thisData[[NUCLEUS_SURFACE_AREA_MIC]][i];
  }
  for(i in 1:length(counts)){
    if(counts[i] > 0){
      averageEmbryoMorphVersusDistance[[colHeading]][i] <- averageEmbryoMorphVersusDistance[[colHeading]][i] + currentEmbryoMorphVersusDistance[i] / counts[i];
    }
  }
}

write.csv(averageEmbryoMorphVersusDistance, paste("outputs", "cell_morphology_versus_distance.csv", sep="/"));