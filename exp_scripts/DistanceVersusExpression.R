source("exp_scripts/Definitions.R");

embryoData <- allData[[embryoHeading]];
embryos <- unique(embryoData);

treatments <- c(CONTROL, TREATED);

averageEmbryoExpressionVersusDistance <- rep(list(rep(list(rep(list(list()), 11)), 2)), 2);
targetCol1 <- paste(NUCLEUS, GATA, MEAN_INTENSITY, sep ="_" );
targetCol2 <- paste(NUCLEUS, YAP, MEAN_INTENSITY, sep ="_" );
targetCol3 <- paste(CYTOPLASM, YAP, MEAN_INTENSITY, sep ="_" );

for(e in embryos){
  thisData <- subset(allData, embryoData==e);
  n <- nrow(thisData);
  t <- thisData[[treatmentHeading]][1];
  colHeading <- paste(treatments[t], MEAN, CELL_VOLUME_MICRONS_VERSUS_DISTANCE, sep="_");
  currentEmbryoExpressionVersusDistance <- matrix(0.0,nrow=11,ncol=2);
  
  counts <- rep(0, 11);
  for(i in 1:n){
    d <- round(thisData[[DISTANCE_TO_CENTRE]][i] * 10) + 1;
    counts[d] <- counts[d] + 1;
    currentEmbryoExpressionVersusDistance[d,1] <- currentEmbryoExpressionVersusDistance[d,1] + thisData[[targetCol1]][i];
    currentEmbryoExpressionVersusDistance[d,2] <- currentEmbryoExpressionVersusDistance[d,2] + thisData[[targetCol2]][i] / thisData[[targetCol3]][i];
  }
  for(i in 1:length(counts)){
    if(counts[i] > 0){
      for(j in 1:length(averageEmbryoExpressionVersusDistance)){
        averageEmbryoExpressionVersusDistance[[j]][[t]][[i]] <- append(averageEmbryoExpressionVersusDistance[[j]][[t]][[i]], currentEmbryoExpressionVersusDistance[i,j] / counts[i]);
      }
    }
  }
}

meansErr <- getMeansErrors(2, 2, averageEmbryoExpressionVersusDistance);

means <- unlist(meansErr[[1]]);
err <- unlist(meansErr[[2]]);

x <- seq(from=0,to=1,by=0.1);
xLabel <- "Normalised Distance to Embryo Centre";
gataExpressLabel <- paste(GATA, "Expression (AU)");
yapExpressLabel <- paste(NUCLEUS, ":", CYTOPLASM, YAP, "Expression Ratio");
gataDistLabel <- paste(NUCLEUS, GATA, "Expression Versus Distance");
yapDistLabel <- paste(NUCLEUS, ":", CYTOPLASM, YAP, "Expression Versus Distance");

savePlot("plots", "gata3_nuc_expression_versus_distance.pdf", x, means[,1], means[,2], err[,1], err[,2], gataDistLabel, xLabel, gataExpressLabel, c(0,60));
savePlot("plots", "nuc_to_cyto_yap1_expression_ratio_versus_distance.pdf", x, means[,3], means[,4], err[,3], err[,4], "", xLabel, yapExpressLabel, c(0,4));

saveMeanErrors(x, means, err, c(xLabel, gataExpressLabel, yapExpressLabel), "outputs", "expression_versus_distance.csv");

