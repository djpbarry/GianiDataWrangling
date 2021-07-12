source("published_scripts/definitions.R");

#simpleImaris <- FALSE;

gtFiles <- list.files(path = gtDir, pattern = "Image Simulator_", recursive = TRUE, full.names = TRUE, include.dirs = TRUE);

#if(simpleImaris){
  IM_CENTROID_X <- "Location_Center_X";
  IM_CENTROID_Y <- "Location_Center_Y";
  IM_CENTROID_Z <- "Location_Center_Z";
#  cpData <- read.csv('Z:/working/barryd/GIANI_Paper/Sim_Data/imaris_batch/stats/Spots_Position_with_Ground_Truth.csv', encoding="UTF-8");
#} else {
#  IM_CENTROID_X <- "Cell.Position.X";
#  IM_CENTROID_Y <- "Cell.Position.Y";
#  IM_CENTROID_Z <- "Cell.Position.Z";
  cpNucData <- read.csv('E:/Dropbox (The Francis Crick)/Debugging/Giani/Sim_Data/CP_Outputs/3d_monolayer_Nuclei_with_Ground_Truth.csv', encoding="UTF-8");
  cpCellData <- read.csv('E:/Dropbox (The Francis Crick)/Debugging/Giani/Sim_Data/CP_Outputs/3d_monolayer_Cells_with_Ground_Truth.csv', encoding="UTF-8");
#}

  cpNucData$Location_Center_X <- cpNucData$Location_Center_X * 0.2929;
  cpNucData$Location_Center_Y <- cpNucData$Location_Center_Y * 0.2929;
  cpNucData$Location_Center_Z <- cpNucData$Location_Center_Z * 1.339;  
  cpCellData$Location_Center_X <- cpCellData$Location_Center_X * 0.2929;
  cpCellData$Location_Center_Y <- cpCellData$Location_Center_Y * 0.2929;
  cpCellData$Location_Center_Z <- cpCellData$Location_Center_Z * 1.339; 
  
allData <- data.frame();

index <- 0;

for (f in gtFiles){
  snrStartPos <- regexpr(SNR, f) + nchar(SNR);
  #snrEndPos <- snrStartPos + regexpr(.Platform$file.sep, substr(f, snrStartPos, nchar(f))) - 2;
  
  snrEndPos <- regexpr('_ncells', f) - 1;
  snr <- as.numeric(substr(f, snrStartPos, snrEndPos));
  
  #runStartPos <- regexpr(RUN, f) + nchar(RUN);
  #runEndPos <- runStartPos + regexpr(.Platform$file.sep, substr(f, runStartPos, nchar(f))) - 2;
  #run <- as.numeric(substr(f, runStartPos, runEndPos));
  
  run <- 1;
  if(grepl("\\(1\\)", f)){
    run <- 2;
  } else if(grepl("\\(2\\)", f)){
    run <- 3;
  }
  
  csvFiles <- list.files(path = f, pattern = ".csv");
  groundTruthData <- read.csv(file.path(f, csvFiles[1]), encoding="UTF-8");
  #groundTruthData <- read.csv(f);
  
  #matchingImarisDir<-imarisDirs[lapply(lapply(imarisDirs, grep, pattern=paste(paste(IM_FOLDER_1, paste(SNR, format(snr, digits=1, nsmall=1), sep=""), format(run, digits=0), sep = "_"), IM_FOLDER_2, sep=" "), value=FALSE), length)>0];
  
  #if(length(matchingImarisDir) < 1){
  #  next;
  #}
  
  #nucPosData <- read.csv(file.path(matchingImarisDir, paste(paste(IM_NUC_POS_FILE_1, paste(SNR, format(snr, digits=1, nsmall=1), sep=""), format(run, digits=0), sep = "_"), IM_NUC_POS_FILE_2, sep=" ")), skip=3);
  #cellVolData <- read.csv(file.path(matchingImarisDir, paste(paste(IM_CELL_VOL_FILE_1, paste(SNR, format(snr, digits=1, nsmall=1), sep=""), format(run, digits=0), sep = "_"), IM_CELL_VOL_FILE_2, sep=" ")), skip=3);
  
  #nucPosData <- cpNucData[grepl(gsub("\\)", "", gsub("\\(", "", paste(basename(f), '_ Sim', sep=" "))), gsub("\\)", "", gsub("\\(", "", cpNucData$Original_Image_Name))),];
  nucPosData <- cpNucData[grepl(basename(f), cpNucData$Original_Image_Name),];
  cellPosData <- cpCellData[grepl(basename(f), cpCellData$Original_Image_Name),];
  
  nIM <- nrow(nucPosData);
  nGT <- nrow(groundTruthData);
  
  if(nIM < 1){
    next;
  }

  appendage <- data.frame(matrix(data = NaN, nrow = nGT, ncol = 11));
  colnames(appendage) <- c(GT_PROP_VOL, IM_CENTROID_X, IM_CENTROID_Y, IM_CENTROID_Z, IM_CELL_VOL, IM_PROP_VOL, NUC_CENTROID_ERROR, CELL_CENTROID_ERROR, VOL_ERROR, NORM_VOL_ERROR, PROP_VOL_ERROR);
  
  thisData <- cbind(groundTruthData, appendage);
  
  for(i in 1:nIM){
    imIndex <- nucPosData$Ground_Truth_Cell_Index[i] - 1;
    gtIndex <- match(imIndex, thisData$Cell_Index);
    imx <- nucPosData[[IM_CENTROID_X]][i];
    imy <- nucPosData[[IM_CENTROID_Y]][i];
    imz <- nucPosData[[IM_CENTROID_Z]][i];

    gtx <- thisData[[GT_CENTROID_X]][gtIndex];
    gty <- thisData[[GT_CENTROID_Y]][gtIndex];
    gtz <- thisData[[GT_CENTROID_Z]][gtIndex];
    dist <- euclidDistance(c(imx,imy,imz), c(gtx, gty, gtz));

    thisData[[IM_CENTROID_X]][gtIndex] <- imx;
    thisData[[IM_CENTROID_Y]][gtIndex] <- imy;
    thisData[[IM_CENTROID_Z]][gtIndex] <- imz;
    thisData[[NUC_CENTROID_ERROR]][gtIndex] <- dist;
    
    imIndex <- nucPosData$Ground_Truth_Cell_Index[i] - 1;
    gtIndex <- match(imIndex, thisData$Cell_Index);
    imx <- cellPosData[[IM_CENTROID_X]][i];
    imy <- cellPosData[[IM_CENTROID_Y]][i];
    imz <- cellPosData[[IM_CENTROID_Z]][i];
    
    gtx <- thisData[[GTC_CENTROID_X]][gtIndex];
    gty <- thisData[[GTC_CENTROID_Y]][gtIndex];
    gtz <- thisData[[GTC_CENTROID_Z]][gtIndex];
    dist <- euclidDistance(c(imx,imy,imz), c(gtx, gty, gtz));

    thisData[[CELL_CENTROID_ERROR]][gtIndex] <- dist;
  }

  extraCols <- cbind(matrix(data=snr,ncol=1,nrow=nGT),
                     matrix(data=run,ncol=1,nrow=nGT),
                     matrix(data=nGT,ncol=1,nrow=nGT),
                     matrix(data=nIM,ncol=1,nrow=nGT),
                     matrix(data=nIM - nGT,ncol=1,nrow=nGT),
                     matrix(data=(nIM - nGT) / nGT,ncol=1,nrow=nGT),
                     matrix(data=index,ncol=1,nrow=nGT));
  colnames(extraCols) <- c(SNR, RUN, GROUND_TRUTH_N, MEASURED_N, CELL_COUNT_ERROR, PROP_CELL_COUNT_ERROR, INDEX);
  
  thisData <- cbind(thisData, extraCols);
  allData <- rbind(allData, thisData);
  
  index <- index + 1;
}

filename <- "all_cell_profiler_data.csv";

write.csv(allData, file.path("outputs", filename));
