source("sim_scripts/Definitions.R");

nuc <- 1;
cyto <- 3;
cell <- 5;

files <- list.files(path = directory, pattern = "Image Simulator_Output", recursive = TRUE, full.names = TRUE, include.dirs = TRUE);

allData <- data.frame();
falsePositives <- data.frame();

for (f in files){
  groundTruthData <- read.csv(file.path(f, GROUND_TRUTH));
  gianiData <- read.csv(file.path(f, GIANI));
  gianiColNames <- colnames(gianiData);
  nucGianiData <- subset(gianiData, grepl(NUCLEUS, gianiData[[LABEL]]));
  cellGianiData <- subset(gianiData, grepl(CELL, gianiData[[LABEL]]));
  nGD <- nrow(nucGianiData);
  nGT <- nrow(groundTruthData);
  gtFoundCol <- matrix(data=0,ncol=1,nrow=nGD);
  colnames(gtFoundCol) <- c(GROUND_TRUTH_FOUND);
  nucGianiData <- cbind(nucGianiData, gtFoundCol);
  appendage <- data.frame(matrix(data = NaN, nrow = nGT, ncol = 5))
  colnames(appendage) <- c(GD_CENTROID_X, GD_CENTROID_Y, GD_CENTROID_Z, CENTROID_ERROR,gianiColNames[5]);
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
    
    if(!is.nan(thisData[[CENTROID_ERROR]][minDistIndex])){
      print("duplicate");
    }
    if(nucGianiData[[GROUND_TRUTH_FOUND]][i] < 1 && (is.nan(thisData[[CENTROID_ERROR]][minDistIndex]) || thisData[[CENTROID_ERROR]][minDistIndex] > minDist)){
      thisData[[GD_CENTROID_X]][minDistIndex] <- gdx;
      thisData[[GD_CENTROID_Y]][minDistIndex] <- gdy;
      thisData[[GD_CENTROID_Z]][minDistIndex] <- gdz;
      thisData[[CENTROID_ERROR]][minDistIndex] <- minDist;
      thisData[[gianiColNames[5]]][minDistIndex] <- cellGianiData[[gianiColNames[5]]][i];
      nucGianiData[[GROUND_TRUTH_FOUND]][i] <- 1;
    } else {
      #thisFalsePositive <- data.frame(c(gdx, gdy, gdz, cellGianiData[[gianiColNames[5]]][i]));
      falsePositives <- rbind(falsePositives, data.frame(c(gdx, gdy, gdz, cellGianiData[[gianiColNames[5]]][i])));
    }
  }
  snrStartPos <- regexpr(SNR, f) + nchar(SNR);
  snrEndPos <- snrStartPos + regexpr(.Platform$file.sep, substr(f, snrStartPos, nchar(f))) - 2;
  snr <- as.numeric(substr(f, snrStartPos, snrEndPos));
  
  runStartPos <- regexpr(RUN, f) + nchar(RUN);
  runEndPos <- runStartPos + regexpr(.Platform$file.sep, substr(f, runStartPos, nchar(f))) - 2;
  run <- as.numeric(substr(f, runStartPos, runEndPos));
  
  extraCols <- cbind(matrix(data=snr,ncol=1,nrow=nGT), matrix(data=run,ncol=1,nrow=nGT));
  colnames(extraCols) <- c(SNR, RUN);
  thisData <- cbind(thisData, extraCols);
  allData <- rbind(allData, thisData);
}
