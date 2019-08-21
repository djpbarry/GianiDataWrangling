source("imaris_scripts/Definitions.R");

gtDirs <- list.files(path = directory, pattern = "Image Simulator_Output_", recursive = TRUE, full.names = TRUE, include.dirs = TRUE);
imarisDirs <- list.dirs(path = 'Z:/working/barryd/hpc', recursive = TRUE, full.names = TRUE);

allData <- data.frame();

index <- 0;

for (f in gtDirs){
  snrStartPos <- regexpr(SNR, f) + nchar(SNR);
  snrEndPos <- snrStartPos + regexpr(.Platform$file.sep, substr(f, snrStartPos, nchar(f))) - 2;
  snr <- as.numeric(substr(f, snrStartPos, snrEndPos));
  
  runStartPos <- regexpr(RUN, f) + nchar(RUN);
  runEndPos <- runStartPos + regexpr(.Platform$file.sep, substr(f, runStartPos, nchar(f))) - 2;
  run <- as.numeric(substr(f, runStartPos, runEndPos));
  
  groundTruthData <- read.csv(file.path(f, GROUND_TRUTH_FILE));
  
  matchingImarisDir<-imarisDirs[lapply(lapply(imarisDirs, grep, pattern=paste(paste(IM_FOLDER_1, paste(SNR, format(snr, digits=1, nsmall=1), sep=""), format(run, digits=0), sep = "_"), IM_FOLDER_2, sep=" "), value=FALSE), length)>0];
  
  nucPosData <- read.csv(file.path(matchingImarisDir, paste(paste(IM_NUC_POS_FILE_1, paste(SNR, format(snr, digits=1, nsmall=1), sep=""), format(run, digits=0), sep = "_"), IM_NUC_POS_FILE_2, sep=" ")), skip=3);
  cellVolData <- read.csv(file.path(matchingImarisDir, paste(paste(IM_CELL_VOL_FILE_1, paste(SNR, format(snr, digits=1, nsmall=1), sep=""), format(run, digits=0), sep = "_"), IM_CELL_VOL_FILE_2, sep=" ")), skip=3);
  
  nIM <- nrow(nucPosData);
  nGT <- nrow(groundTruthData);

  appendage <- data.frame(matrix(data = NaN, nrow = nGT, ncol = 10))
  colnames(appendage) <- c(GT_PROP_VOL, IM_CENTROID_X, IM_CENTROID_Y, IM_CENTROID_Z, IM_CELL_VOL, IM_PROP_VOL, CENTROID_ERROR, VOL_ERROR, NORM_VOL_ERROR, PROP_VOL_ERROR);
  
  thisData <- cbind(groundTruthData, appendage);
  
  for(i in 1:nIM){  
    imx <- nucPosData[[IM_CENTROID_X]][i];
    imy <- nucPosData[[IM_CENTROID_Y]][i];
    imz <- nucPosData[[IM_CENTROID_Z]][i];
    v1 <- c(imx, imy, imz);
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
      thisData[[IM_CENTROID_X]][minDistIndex] <- imx;
      thisData[[IM_CENTROID_Y]][minDistIndex] <- imy;
      thisData[[IM_CENTROID_Z]][minDistIndex] <- imz;
      thisData[[CENTROID_ERROR]][minDistIndex] <- minDist;
      thisData[[IM_CELL_VOL]][minDistIndex] <- cellVolData[[IM_CELL_VOL]][i];
      thisData[[VOL_ERROR]][minDistIndex] <- thisData[[IM_CELL_VOL]][minDistIndex] - thisData[[GT_CELL_VOLUME_MIC]][minDistIndex];
      thisData[[NORM_VOL_ERROR]][minDistIndex] <- thisData[[VOL_ERROR]][minDistIndex] / thisData[[GT_CELL_VOLUME_MIC]][minDistIndex];
    }
  }
  
  thisData <- subset(thisData, !is.nan(thisData$Cell.Volume));
  
  sumGTVol <- sum(thisData[[GT_CELL_VOLUME_MIC]]);
  sumIMVol <- sum(thisData[[IM_CELL_VOL]]);
  
  for(i in 1:nIM){  
      thisData[[GT_PROP_VOL]][i] <- thisData[[GT_CELL_VOLUME_MIC]][i] / sumGTVol;
      thisData[[IM_PROP_VOL]][i] <- thisData[[IM_CELL_VOL]][i] / sumIMVol;
      thisData[[PROP_VOL_ERROR]][i] <- thisData[[IM_PROP_VOL]][i] - thisData[[GT_PROP_VOL]][i];
  }

  extraCols <- cbind(matrix(data=snr,ncol=1,nrow=nIM),
                     matrix(data=run,ncol=1,nrow=nIM),
                     matrix(data=nGT,ncol=1,nrow=nIM),
                     matrix(data=nIM,ncol=1,nrow=nIM),
                     matrix(data=nIM - nGT,ncol=1,nrow=nIM),
                     matrix(data=index,ncol=1,nrow=nIM));
  colnames(extraCols) <- c(SNR, RUN, GROUND_TRUTH_N, MEASURED_N, CELL_COUNT_ERROR, INDEX);
  
  thisData <- cbind(thisData, extraCols);
  allData <- rbind(allData, thisData);
  
  index <- index + 1;
}

write.csv(allData, paste("outputs", "all_data.csv", sep="/"));