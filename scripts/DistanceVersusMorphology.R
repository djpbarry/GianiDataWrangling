source("scripts/Definitions.R");

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
      for(j in 1:3){
        averageEmbryoMorphVersusDistance[[j]][[t]][[i]] <- append(averageEmbryoMorphVersusDistance[[j]][[t]][[i]], currentEmbryoMorphVersusDistance[i,j] / counts[i]);
      }
    }
  }
}

# Nuclear volume plot

means <- matrix(nrow=11,ncol=6);
err <- matrix(nrow=11,ncol=6);

for(j in 1:3){
  for(t in 1:2){
    for(i in 1:length(averageEmbryoMorphVersusDistance[[j]][[t]])){
      if(length(averageEmbryoMorphVersusDistance[[j]][[t]][[i]]) > 1){
        means[i, t + (j-1)*2] <- mean(unlist(averageEmbryoMorphVersusDistance[[j]][[t]][[i]]));
        interval <- CI(unlist(averageEmbryoMorphVersusDistance[[j]][[t]][[i]]), confLevel);
        err[i,t+(j-1)*2] <- interval[1] - interval[2];
      } else {
        means[i,t+(j-1)*2] <- unlist(averageEmbryoMorphVersusDistance[[j]][[t]][[i]]);
        err[i,t+(j-1)*2] <- 0.0;
      }
    }
  }
}

x <- seq(from=0,to=1,by=0.1);

pdf(paste("plots", "nuc_volume_versus_distance.pdf", sep=.Platform$file.sep));

plot(x, means[,1], xlab="Normalised Distance from Embryo Centre", main="Nuclear Volume Versus Distance", ylab=expression(paste("Nuclear Volume (",mu,"m^3)")), ylim=c(500,3000), col="red", pch=15);
arrows(x, means[,1]-err[,1], x, means[,1]+err[,1], length=0.05, angle=90, code=3, col="red");

points(x, means[,2], col="green", pch=15);
arrows(x, means[,2]-err[,2], x, means[,2]+err[,2], length=0.05, angle=90, code=3, col="green");

legend("bottomright", c(CONTROL, TREATED), fill=c("red", "green"));

dev.off();


#Cell Volume Plot

pdf(paste("plots", "cell_volume_versus_distance.pdf", sep=.Platform$file.sep));

plot(x, means[,3], xlab="Normalised Distance from Embryo Centre", main="Cell Volume Versus Distance", ylab=expression(paste("Cell Volume (",mu,"m^3)")),  ylim=c(500,20000), col="red", pch=15);
arrows(x, means[,3]-err[,3], x, means[,3]+err[,3], length=0.05, angle=90, code=3, col="red");

points(x, means[,4], col="green", pch=15);
arrows(x, means[,4]-err[,4], x, means[,4]+err[,4], length=0.05, angle=90, code=3, col="green");

legend("bottomright", c(CONTROL, TREATED), fill=c("red", "green"));

dev.off();

#Nuclear:Cell Volume Ratio Plot

pdf(paste("plots", "nuc_to_cell_volume_ratio_versus_distance.pdf", sep=.Platform$file.sep));

plot(x, means[,5], xlab="Normalised Distance from Embryo Centre", main="Nuclear:Cell Volume Ratio Versus Distance", ylab="Nuclear:Cell Volume Ratio",  ylim=c(0,0.5), col="red", pch=15);
arrows(x, means[,5]-err[,5], x, means[,5]+err[,5], length=0.05, angle=90, code=3, col="red");

points(x, means[,6], col="green", pch=15);
arrows(x, means[,6]-err[,6], x, means[,6]+err[,6], length=0.05, angle=90, code=3, col="green");

legend("bottomright", c(CONTROL, TREATED), fill=c("red", "green"));

dev.off();

# write.csv(averageEmbryoMorphVersusDistance, paste("outputs", "cell_morphology_versus_distance.csv", sep="/"));