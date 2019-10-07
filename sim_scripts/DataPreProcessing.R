source("sim_scripts/Definitions.R");

gtFiles <- list.files(path = gtDir, pattern = ".csv", recursive = TRUE, full.names = TRUE);
gianiFiles <- list.files(path = gianiDir, pattern = ".csv", recursive = TRUE, full.names = TRUE);

allData <- data.frame();

index <- 0;

for (fileIndex in 1:length(gtFiles)){
  groundTruthData <- read.csv(gtFiles[fileIndex]);
  gianiData <- read.csv(gianiFiles[fileIndex]);
  
  gianiColNames <- colnames(gianiData);
  nucGianiData <- subset(gianiData, grepl(NUCLEUS, gianiData[[LABEL]]));
  cellGianiData <- subset(gianiData, grepl(CELL, gianiData[[LABEL]]));
  
  nGD <- nrow(nucGianiData);
  nGT <- nrow(groundTruthData);
  
  gtFoundCol <- matrix(data=0,ncol=1,nrow=nGD);
  colnames(gtFoundCol) <- c(GROUND_TRUTH_FOUND);
  nucGianiData <- cbind(nucGianiData, gtFoundCol);
  
  appendage <- data.frame(matrix(data = NaN, nrow = nGT, ncol = 10))
  colnames(appendage) <- c(GT_PROP_VOL, GD_CENTROID_X, GD_CENTROID_Y, GD_CENTROID_Z,  gianiColNames[5], GD_PROP_VOL, CENTROID_ERROR, VOL_ERROR, NORM_VOL_ERROR, PROP_VOL_ERROR);
  
  thisData <- cbind(groundTruthData, appendage);
  
  for(i in 1:nGD){  
    gdx <- gianiData[[GD_CENTROID_X]][i];
    gdy <- gianiData[[GD_CENTROID_Y]][i];
    gdz <- gianiData[[GD_CENTROID_Z]][i];
    v1 <- c(gdx, gdy, gdz);
    minDist <- .Machine$double.xmax;
    minDistIndex <- -1;
    for(j in 1:nGT){
      gtx <- thisData[[GT_CENTROID_X]][j];
      gty <- thisData[[GT_CENTROID_Y]][j];
      gtz <- thisData[[GT_CENTROID_Z]][j];
      dist <- euclidDistance(v1, c(gtx, gty, gtz));
      if(dist < minDist){
        minDist <- dist;
        minDistIndex <- j;
      }
    }
    if(is.nan(thisData[[CENTROID_ERROR]][minDistIndex]) || thisData[[CENTROID_ERROR]][minDistIndex] > minDist){
      thisData[[GD_CENTROID_X]][minDistIndex] <- gdx;
      thisData[[GD_CENTROID_Y]][minDistIndex] <- gdy;
      thisData[[GD_CENTROID_Z]][minDistIndex] <- gdz;
      thisData[[CENTROID_ERROR]][minDistIndex] <- minDist;
      thisData[[gianiColNames[5]]][minDistIndex] <- cellGianiData[[gianiColNames[5]]][i];
      thisData[[VOL_ERROR]][minDistIndex] <- abs(thisData[[gianiColNames[5]]][minDistIndex] - thisData[[GT_CELL_VOLUME_MIC]][minDistIndex]);
      thisData[[NORM_VOL_ERROR]][minDistIndex] <- thisData[[VOL_ERROR]][minDistIndex] / thisData[[GT_CELL_VOLUME_MIC]][minDistIndex];
      nucGianiData[[GROUND_TRUTH_FOUND]][i] <- 1;
    }
  }
  
  sumGTVol <- sum(thisData[[GT_CELL_VOLUME_MIC]]);
  sumGDVol <- sum(thisData[[gianiColNames[5]]]);
  
  for(i in 1:nGT){  
      thisData[[GT_PROP_VOL]][i] <- thisData[[GT_CELL_VOLUME_MIC]][i] / sumGTVol;
      thisData[[GD_PROP_VOL]][i] <- thisData[[gianiColNames[5]]][i] / sumGDVol;
      thisData[[PROP_VOL_ERROR]][i] <- thisData[[GD_PROP_VOL]][i] - thisData[[GT_PROP_VOL]][i];
  }
  
  snrStartPos <- regexpr(SNR, gtFiles[fileIndex]) + nchar(SNR);
  snrEndPos <- snrStartPos + regexpr("_", substr(gtFiles[fileIndex], snrStartPos, nchar(gtFiles[fileIndex]))) - 2;
  snr <- as.numeric(substr(gtFiles[fileIndex], snrStartPos, snrEndPos));
  
  extraCols <- cbind(matrix(data=snr,ncol=1,nrow=nGT),
                     matrix(data=nGT,ncol=1,nrow=nGT),
                     matrix(data=nGD,ncol=1,nrow=nGT),
                     matrix(data=nGD - nGT,ncol=1,nrow=nGT),
                     matrix(data=(nGD - nGT)/nGT,ncol=1,nrow=nGT),
                     matrix(data=index,ncol=1,nrow=nGT));
  colnames(extraCols) <- c(SNR, GROUND_TRUTH_N, MEASURED_N, CELL_COUNT_ERROR, PROP_CELL_COUNT_ERROR, INDEX);
  
  thisData <- cbind(thisData, extraCols);
  allData <- rbind(allData, thisData);
  
  index <- index + 1;
}

write.csv(allData, paste("outputs", "sim_all_data.csv", sep="/"));