source("exp_scripts/Definitions.R");

embryoData <- allData[[embryoHeading]];
embryos <- unique(embryoData);

treatments <- c(CONTROL, TREATED);

averageEmbryoMorphVersusDistance <- rep(list(rep(list(rep(list(list()), 11)), 2)), 3);

for(e in embryos){
  thisData <- subset(allData, embryoData==e);
  n <- nrow(thisData);
  t <- thisData[[treatmentHeading]][1];
  colHeading <- paste(treatments[t], MEAN, CELL_VOLUME_MICRONS_VERSUS_DISTANCE, sep="_");
  currentEmbryoMorphVersusDistance <- matrix(0.0,nrow=11,ncol=3);
  
  counts <- rep(0, 11);
  for(i in 1:n){
    d <- round(thisData[[DISTANCE_TO_CENTRE]][i] * 10) + 1;
    counts[d] <- counts[d] + 1;
    currentEmbryoMorphVersusDistance[d,1] <- currentEmbryoMorphVersusDistance[d,1] + thisData[[NUCLEUS_VOLUME_MIC]][i];
    currentEmbryoMorphVersusDistance[d,2] <- currentEmbryoMorphVersusDistance[d,2] + thisData[[CELL_VOLUME_MIC]][i];
    currentEmbryoMorphVersusDistance[d,3] <- currentEmbryoMorphVersusDistance[d,3] + thisData[[NUCLEUS_VOLUME_MIC]][i] / thisData[[CELL_VOLUME_MIC]][i];
  }
  for(i in 1:length(counts)){
    if(counts[i] > 0){
      for(j in 1:length(averageEmbryoMorphVersusDistance)){
        averageEmbryoMorphVersusDistance[[j]][[t]][[i]] <- append(averageEmbryoMorphVersusDistance[[j]][[t]][[i]], currentEmbryoMorphVersusDistance[i,j] / counts[i]);
      }
    }
  }
}

meansErr <- getMeansErrors(3, 2, averageEmbryoMorphVersusDistance);

means <- unlist(meansErr[[1]]);
err <- unlist(meansErr[[2]]);

x <- seq(from=0,to=1,by=0.1);
xLabel <- "Normalised Distance to Embryo Centre";
nucVolLabel <- paste(NUCLEUS, "Volume (microns^3)");
nucVolHeading <- paste(NUCLEUS, "Volume Versus Distance");
cellVolLabel <- paste(CELL, "Volume (microns^3)");
cellVolHeading <- paste(CELL, "Volume Versus Distance");
nucCellRatioLabel <- paste(CELL, "Volume (microns^3)");
nucCellRatioHeading <- paste(CELL, "Volume Versus Distance");

savePlot("plots", "nuc_volume_versus_distance.pdf", x, means[,1], means[,2], err[,1], err[,2], nucVolHeading, xLabel, nucVolLabel, c(500,3000));
savePlot("plots", "cell_volume_versus_distance.pdf", x, means[,3], means[,4], err[,3], err[,4], cellVolHeading, xLabel, cellVolLabel, c(500,20000));
savePlot("plots", "nuc_to_cell_volume_ratio_versus_distance.pdf", x, means[,5], means[,6], err[,5], err[,6], nucCellRatioHeading, xLabel, nucCellRatioLabel, c(0,0.5));

saveMeanErrors(x, means, err, c(xLabel, nucVolLabel,cellVolLabel,nucCellRatioLabel), "outputs", "cell_morphology_versus_distance.csv");
