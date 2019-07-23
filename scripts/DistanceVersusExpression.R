source("scripts/Definitions.R");

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
      for(j in 1:2){
        averageEmbryoExpressionVersusDistance[[j]][[t]][[i]] <- append(averageEmbryoExpressionVersusDistance[[j]][[t]][[i]], currentEmbryoExpressionVersusDistance[i,j] / counts[i]);
      }
    }
  }
}

# Nuclear volume plot

means <- matrix(nrow=11,ncol=4);
err <- matrix(nrow=11,ncol=4);

for(j in 1:2){
  for(t in 1:2){
    for(i in 1:length(averageEmbryoExpressionVersusDistance[[j]][[t]])){
      if(length(averageEmbryoExpressionVersusDistance[[j]][[t]][[i]]) > 1){
        means[i, t + (j-1)*2] <- mean(unlist(averageEmbryoExpressionVersusDistance[[j]][[t]][[i]]));
        interval <- CI(unlist(averageEmbryoExpressionVersusDistance[[j]][[t]][[i]]), confLevel);
        err[i,t+(j-1)*2] <- interval[1] - interval[2];
      } else {
        means[i,t+(j-1)*2] <- unlist(averageEmbryoExpressionVersusDistance[[j]][[t]][[i]]);
        err[i,t+(j-1)*2] <- 0.0;
      }
    }
  }
}

x <- seq(from=0,to=1,by=0.1);

pdf(paste("plots", "gata3_nuc_expression_versus_distance.pdf", sep=.Platform$file.sep));

plot(x, means[,1], xlab="Normalised Distance from Embryo Centre", main="Nuclear GATA3 Expression Versus Distance", ylab="GATA3 Expression (AU)", ylim=c(0,60), col="red", pch=15);
arrows(x, means[,1]-err[,1], x, means[,1]+err[,1], length=0.05, angle=90, code=3, col="red");

points(x, means[,2], col="green", pch=15);
arrows(x, means[,2]-err[,2], x, means[,2]+err[,2], length=0.05, angle=90, code=3, col="green");

legend("topright", c(CONTROL, TREATED), fill=c("red", "green"));

dev.off();



pdf(paste("plots", "nuc_to_cyto_yap1_expression_ratio_versus_distance.pdf", sep=.Platform$file.sep));

plot(x, means[,3], xlab="Normalised Distance from Embryo Centre", main="Nuclear:Cytoplasm YAP1 Expression Versus Distance", ylab="Nuclear:Cytoplasm YAP1 Expression Ratio", ylim=c(0,4), col="red", pch=15);
arrows(x, means[,3]-err[,3], x, means[,3]+err[,3], length=0.05, angle=90, code=3, col="red");

points(x, means[,4], col="green", pch=15);
arrows(x, means[,4]-err[,4], x, means[,4]+err[,4], length=0.05, angle=90, code=3, col="green");

legend("bottomright", c(CONTROL, TREATED), fill=c("red", "green"));

dev.off();


# write.csv(averageEmbryoExpressionVersusDistance, paste("outputs", "cell_morphology_versus_distance.csv", sep="/"));