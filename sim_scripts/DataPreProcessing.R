source("sim_scripts/Definitions.R");


gtFiles <-
  sort(list.files(
    path = gtDir,
    pattern = ".csv",
    recursive = TRUE,
    full.names = TRUE
  ))

#gianiDirs <- list.files(path = gianiDir, pattern = "GIANI v2.068_Sim_Image_snr[[:digit:]]{1}.[[:digit:]]{6}_ncells[[:digit:]]{2}.tif_S0_Output_[[:print:]]", recursive = TRUE, full.names = TRUE, include.dirs = TRUE);
gianiFiles <-
  sort(list.files(
    path = gianiDir,
    pattern = "_with_Ground_Truth.csv",
    recursive = TRUE,
    full.names = TRUE
  ))

allData <- data.frame()

index <- 0


for (fileIndex in 1:length(gtFiles)) {
  groundTruthData <- read.csv(gtFiles[fileIndex])
  groundTruthData1 <- groundTruthData;
  #correct ordering cell centroids in ground truth data
  
  groundTruthData$Cell_Centroid_X <- groundTruthData$Cell_Centroid_X[groundTruthData$Cell_Index + 1];
  groundTruthData$Cell_Centroid_Y <- groundTruthData$Cell_Centroid_Y[groundTruthData$Cell_Index + 1];
  groundTruthData$Cell_Centroid_Z <- groundTruthData$Cell_Centroid_Z[groundTruthData$Cell_Index + 1];
  
  run <- 1;
  if(grepl("\\(1\\)", gtFiles[fileIndex])){
    run <- 2;
  } else if(grepl("\\(2\\)", gtFiles[fileIndex])){
    run <- 3;
  }
  
  #print(gtFiles[fileIndex]);
  gianiData <- read.csv(gianiFiles[2 * fileIndex - 1])
  
  gianiColNames <- colnames(gianiData)
  
  nucGianiData <-
    subset(gianiData, grepl(NUCLEUS, gianiData[[LABEL]]))
  
  cellGianiData <-
    subset(gianiData, grepl(CELL, gianiData[[LABEL]]))
  
  nGD <- nrow(nucGianiData)
  
  nGT <- nrow(groundTruthData)
  
  gtFoundCol <- matrix(data = 0,
                       ncol = 1,
                       nrow = nGD)
  
  colnames(gtFoundCol) <- c(GROUND_TRUTH_FOUND)
  
  nucGianiData <- cbind(nucGianiData, gtFoundCol)
  
  appendage <- data.frame(matrix(
    data = NaN,
    nrow = nGT,
    ncol = 15
  ))
  colnames(appendage) <-
    c(
      GT_PROP_VOL,
      M_NUC_CENTROID_X,
      M_NUC_CENTROID_Y,
      M_NUC_CENTROID_Z,
      M_CELL_CENTROID_X,
      M_CELL_CENTROID_Y,
      M_CELL_CENTROID_Z,
      gianiColNames[5],
      GD_PROP_VOL,
      NUC_CENTROID_ERROR,
      CELL_CENTROID_ERROR,
      VOL_ERROR,
      NORM_VOL_ERROR,
      PROP_VOL_ERROR,
      NORM_DIST
    )
  
  thisData <- cbind(groundTruthData, appendage)
  
  for (i in 1:nGD) {
    gianiIndex <- nucGianiData$Ground_Truth_Cell_Index[i] - 1;
    gtIndex <- match(gianiIndex, thisData$Cell_Index);  
    
    gdx <- gianiData[[GD_CENTROID_X]][i]
    
    gdy <- gianiData[[GD_CENTROID_Y]][i]
    
    gdz <- gianiData[[GD_CENTROID_Z]][i]
    
    gdcx <- cellGianiData[[GD_CENTROID_X]][i]
    
    gdcy <- cellGianiData[[GD_CENTROID_Y]][i]
    
    gdcz <- cellGianiData[[GD_CENTROID_Z]][i]
    
    #v1 <- c(gdx, gdy, gdz)
    
    #minDist <- .Machine$double.xmax
    
    #minDistIndex <- -1
    
    #for (j in 1:nGT) {
      
      gtx <- thisData[[GT_CENTROID_X]][gtIndex]
      
      gty <- thisData[[GT_CENTROID_Y]][gtIndex]
      
      gtz <- thisData[[GT_CENTROID_Z]][gtIndex]
      
      dist <- euclidDistance(c(gdx, gdy, gdz), c(gtx, gty, gtz))
      
     # if (dist < minDist) {
      #  minDist <- dist
        
       # minDistIndex <- j
        
    #  }
    #}
   # if (is.nan(thisData[[NUC_CENTROID_ERROR]][gtIndex]) ||
    #    thisData[[NUC_CENTROID_ERROR]][gtIndex] > minDist) {
      thisData[[M_NUC_CENTROID_X]][gtIndex] <- gdx
      
      thisData[[M_NUC_CENTROID_Y]][gtIndex] <- gdy
      
      thisData[[M_NUC_CENTROID_Z]][gtIndex] <- gdz
      
      thisData[[NUC_CENTROID_ERROR]][gtIndex] <- dist
      
      thisData[[M_CELL_CENTROID_X]][gtIndex] <- gdcx
      
      thisData[[M_CELL_CENTROID_Y]][gtIndex] <- gdcy
      
      thisData[[M_CELL_CENTROID_Z]][gtIndex] <- gdcz
      
      thisData[[CELL_CENTROID_ERROR]][gtIndex] <-
        euclidDistance(
          c(gdcx, gdcy, gdcz),
          c(thisData[[GTC_CENTROID_X]][gtIndex],
            thisData[[GTC_CENTROID_Y]][gtIndex],
            thisData[[GTC_CENTROID_Z]][gtIndex])
        )
      
      thisData[[gianiColNames[5]]][gtIndex] <-
        cellGianiData[[gianiColNames[5]]][i]
      
      thisData[[VOL_ERROR]][gtIndex] <-
        thisData[[gianiColNames[5]]][gtIndex] - thisData[[GT_CELL_VOLUME_MIC]][gtIndex]
      
      thisData[[NORM_VOL_ERROR]][gtIndex] <-
        thisData[[VOL_ERROR]][gtIndex] / thisData[[GT_CELL_VOLUME_MIC]][gtIndex]
      
      nucGianiData[[GROUND_TRUTH_FOUND]][i] <- 1
      
      thisData[[NORM_DIST]] <- gianiData[[NORM_DIST]][i]
      
    #}
  }
  
  sumGTVol <-
    sum(thisData[[GT_CELL_VOLUME_MIC]][!is.nan(thisData[[GT_CELL_VOLUME_MIC]])])
  
  sumGDVol <-
    sum(thisData[[gianiColNames[5]]][!is.nan(thisData[[gianiColNames[5]]])])
  
  normFactor <- sumGTVol / sumGDVol
  
  for (i in 1:nGT) {
    thisData[[GT_PROP_VOL]][i] <-
      thisData[[GT_CELL_VOLUME_MIC]][i] / sumGTVol
    
    thisData[[GD_PROP_VOL]][i] <-
      thisData[[gianiColNames[5]]][i] / sumGDVol
    
    thisData[[PROP_VOL_ERROR]][i] <-
      normFactor * thisData[[GD_PROP_VOL]][i] - thisData[[GT_PROP_VOL]][i]
    
  }
  
  snrStartPos <- regexpr(SNR, gtFiles[fileIndex]) + nchar(SNR)
  
  snrEndPos <-
    snrStartPos + regexpr("_", substr(gtFiles[fileIndex], snrStartPos, nchar(gtFiles[fileIndex]))) - 2
  
  snr <-
    as.numeric(substr(gtFiles[fileIndex], snrStartPos, snrEndPos))
  
  extraCols <- cbind(
    matrix(
      data = snr,
      ncol = 1,
      nrow = nGT
    ),
    matrix(data=run,ncol=1,nrow=nGT),
    matrix(
      data = nGT,
      ncol = 1,
      nrow = nGT
    ),
    matrix(
      data = nGD,
      ncol = 1,
      nrow = nGT
    ),
    matrix(
      data = nGD - nGT,
      ncol = 1,
      nrow = nGT
    ),
    matrix(
      data = (nGD - nGT) / nGT,
      ncol = 1,
      nrow = nGT
    ),
    matrix(
      data = index,
      ncol = 1,
      nrow = nGT
    )
  )
  
  colnames(extraCols) <-
    c(SNR,
      RUN,
      GROUND_TRUTH_N,
      MEASURED_N,
      CELL_COUNT_ERROR,
      PROP_CELL_COUNT_ERROR,
      INDEX)
  
  
  thisData <- cbind(thisData, extraCols)
  
  allData <- rbind(allData, thisData)
  
  index <- index + 1
  
}

write.csv(allData, "E:/Dropbox (The Francis Crick)/Giani_Imaris_Comp/CompiledGianiData.csv");
