source("published_scripts/definitions.R");

genPlot <- function(cellCountErrors, x, y, z, midrange, rangelimits, rangebreaks, rangename, labels){
  red <- "#DD0000";
  green <- "#008800";
  yellow <- "#FFFF00";
  blue <- "#0000FF";
  grey <- "#999999";
  orange <- "#FF6600";
  
  #pdf(paste("plots", "sim_cell_count_errors.pdf", sep=.Platform$file.sep));
  
  #dev.off();
  axislabel <- element_text(hjust=0.5, size=25, colour = "black");
  
  p <- ggplot(cellCountErrors, aes(as.character(.data[[x]]), .data[[y]]));
  p <- p + geom_raster(aes(fill = .data[[z]]),interpolate=FALSE);
  p <- p + scale_fill_gradient2(low=blue, high=orange, mid = grey, midpoint = midrange, name = rangename, limits = rangelimits, breaks =rangebreaks);
  p <- p + xlab(labels[1]) + ylab(labels[2]) + ggtitle(labels[3]);
  p <- p + scale_x_discrete(labels=c("0.10", "0.25", "0.50", "0.75", "1.00", "2.00", "3.00"));
  p <- p + theme_minimal();
  p + theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel, plot.title = axislabel);
  
}

genErrorData <- function(allData, imaris){
  snrs <- unique(allData[[SNR]]);
  cellCounts <- unique(allData[[GROUND_TRUTH_N]]);
  
  snrs<-sort(snrs);
  cellCounts<-sort(cellCounts);
  
  cellCountErrors <- data.frame();
  
  for(i in snrs){
    temp1<-allData[allData[[SNR]] == i, ];
    for(j in cellCounts){
      temp2 <- temp1[temp1$Ground_Truth_Cell_Count == j, ];
      runs <- unique(temp2$run_);
      errors <- matrix(0,3,4);
      for(r in runs){
        temp3 <- temp2[temp2$run_ == r, ];
        #temp4 <- temp3[is.finite(temp3$Centroid_Error),];
        temp4 <- temp3[is.finite(temp3$Nuc_Centroid_Error),];
        errors[r,1] <- mean(temp4$Cell_Count_Error);
        errors[r,2] <- mean(temp4$Proportional_Cell_Count_Error);
        #errors[r,3] <- mean(temp4$Centroid_Error);
        errors[r,3] <- mean(temp4$Nuc_Centroid_Error);
        if(!imaris){
          errors[r,4] <- mean(temp4$Cell_Centroid_Error);
        }
      }
      
      if(imaris) {
        cellCountErrors <- rbind(cellCountErrors, c(i, j, mean(abs(errors[,1])), mean(abs(errors[,2])), mean(errors[,3])));
      } else {
        cellCountErrors <- rbind(cellCountErrors, c(i, j, mean(abs(errors[,1])), mean(abs(errors[,2])), mean(errors[,3]), mean(errors[,4])));
      }
    } 
  }
  if(imaris){
    colnames(cellCountErrors) <- c(SNR, N_CELLS, CELL_COUNT_ERROR, "PROP_CELL_COUNT_ERROR", MEAN_NUC_CENTROID_ERROR);
  } else {
    colnames(cellCountErrors) <- c(SNR, N_CELLS, CELL_COUNT_ERROR, "PROP_CELL_COUNT_ERROR", MEAN_NUC_CENTROID_ERROR, MEAN_CELL_CENTROID_ERROR);
  }
  return(cellCountErrors);
}

simpleGianiData <- read.csv(file.path("outputs", "all_simple_giani_data.csv"));
advancedGianiData <- read.csv(file.path("outputs", "all_advanced_giani_data.csv"));
simpleImarisData <- read.csv(file.path("outputs", "all_simple_imaris_data.csv"));
advancedImarisData <- read.csv(file.path("outputs", "all_simple_imaris_data.csv"));

simpleGianiErrors <- genErrorData(simpleGianiData, FALSE);
advancedGianiErrors <- genErrorData(advancedGianiData, FALSE);
simpleImarisErrors <- genErrorData(simpleImarisData, TRUE);
advancedImarisErrors <- genErrorData(advancedImarisData, TRUE);

#cellCountMin =  floor(min(simpleGianiErrors$Cell_Count_Error, advancedGianiErrors$Cell_Count_Error));
#cellCountMax =  ceiling(max(simpleGianiErrors$Cell_Count_Error, advancedGianiErrors$Cell_Count_Error));
cellCountMin =  round(min(advancedGianiErrors$Cell_Count_Error, advancedImarisErrors$Cell_Count_Error));
cellCountMax =  round(max(advancedGianiErrors$Cell_Count_Error, advancedImarisErrors$Cell_Count_Error));

cellCountMid <- cellCountMin + (cellCountMax - cellCountMin) / 2;
cellCountBreaks <- c(cellCountMin, cellCountMid, cellCountMax);
cellCountLimits <- c(cellCountMin, cellCountMax);

#nucLocMin = floor(10*min(simpleGianiErrors$Mean_Nuclear_Centroid_Error, advancedGianiErrors$Mean_Nuclear_Centroid_Error))/10.0;
#nucLocMax = ceiling(10*max(simpleGianiErrors$Mean_Nuclear_Centroid_Error, advancedGianiErrors$Mean_Nuclear_Centroid_Error))/10.0;
nucLocMin = floor(10*min(advancedGianiErrors$Mean_Nuclear_Centroid_Error, advancedImarisErrors$Mean_Nuclear_Centroid_Error))/10.0;
nucLocMax = ceiling(10*max(advancedGianiErrors$Mean_Nuclear_Centroid_Error, advancedImarisErrors$Mean_Nuclear_Centroid_Error))/10.0;

nucLocMid <- nucLocMin + (nucLocMax - nucLocMin) / 2;
nucLocBreaks <- c(nucLocMin, nucLocMid, nucLocMax);
nucLocLimits <- c(nucLocMin, nucLocMax);

cellLocMin =  floor(min(simpleGianiErrors$Mean_Cell_Centroid_Error, advancedGianiErrors$Mean_Cell_Centroid_Error));
cellLocMax = ceiling(max(simpleGianiErrors$Mean_Cell_Centroid_Error, advancedGianiErrors$Mean_Cell_Centroid_Error));
#cellLocMin =  round(min(advancedGianiErrors$Mean_Cell_Centroid_Error, advancedImarisErrors$Mean_Cell_Centroid_Error));
#cellLocMax =  round(max(advancedGianiErrors$Mean_Cell_Centroid_Error, advancedImarisErrors$Mean_Cell_Centroid_Error));

cellLocMid <- cellLocMin + (cellLocMax - cellLocMin) / 2;
cellLocBreaks <- c(cellLocMin, cellLocMid, cellLocMax);
cellLocLimits <- c(cellLocMin, cellLocMax);

genPlot(simpleGianiErrors, as.character(SNR), N_CELLS, CELL_COUNT_ERROR, cellCountMid, cellCountLimits, cellCountBreaks, expression(paste("E" [c])), c("SNR", "Number of Cells", "GIANI Cell Count Error"));
genPlot(simpleGianiErrors, as.character(SNR), N_CELLS, MEAN_NUC_CENTROID_ERROR, nucLocMid, nucLocLimits, nucLocBreaks, expression(paste("E" [n])), c("SNR", "Number of Cells", "GIANI Nuclear Localisation Error (\U03BCm)"));
genPlot(simpleGianiErrors, as.character(SNR), N_CELLS, MEAN_CELL_CENTROID_ERROR, cellLocMid, cellLocLimits, cellLocBreaks, expression(paste("E" [n])), c("SNR", "Number of Cells", "GIANI Cell Localisation Error (\U03BCm)"));

genPlot(advancedGianiErrors, as.character(SNR), N_CELLS, CELL_COUNT_ERROR, cellCountMid, cellCountLimits, cellCountBreaks, expression(paste("E" [c])), c("SNR", "Number of Cells", "GIANI Cell Count Error"));
genPlot(advancedGianiErrors, as.character(SNR), N_CELLS, MEAN_NUC_CENTROID_ERROR, nucLocMid, nucLocLimits, nucLocBreaks, expression(paste("E" [n])), c("SNR", "Number of Cells", "GIANI Nuclear Localisation Error (\U03BCm)"));
genPlot(advancedGianiErrors, as.character(SNR), N_CELLS, MEAN_CELL_CENTROID_ERROR, cellLocMid, cellLocLimits, cellLocBreaks, expression(paste("E" [n])), c("SNR", "Number of Cells", "GIANI Cell Localisation Error (\U03BCm)"));

#genPlot(simpleImarisErrors, as.character(SNR), N_CELLS, CELL_COUNT_ERROR, cellCountMid, cellCountLimits, cellCountBreaks, expression(paste("E" [c])), c("SNR", "Number of Cells", "Cell Count Error"));
#genPlot(simpleImarisErrors, as.character(SNR), N_CELLS, MEAN_NUC_CENTROID_ERROR, nucLocMid, nucLocLimits, nucLocBreaks, expression(paste("E" [n])), c("SNR", "Number of Cells", "Localisation Error"));

genPlot(advancedImarisErrors, as.character(SNR), N_CELLS, CELL_COUNT_ERROR, cellCountMid, cellCountLimits, cellCountBreaks, expression(paste("E" [c])), c("SNR", "Number of Cells", "Imaris Cell Count Error"));
genPlot(advancedImarisErrors, as.character(SNR), N_CELLS, MEAN_NUC_CENTROID_ERROR, nucLocMid, nucLocLimits, nucLocBreaks, expression(paste("E" [n])), c("SNR", "Number of Cells", "Imaris Nuclear Localisation Error (\U03BCm)"));


#Centroid error plot
p <- ggplot(cellCountErrors, aes(as.character(snr), N_CELLS));
p <- p + geom_raster(aes(fill = MEAN_NUC_CENTROID_ERROR),interpolate=FALSE);
p <- p + scale_fill_gradient2(low=green, high=red, mid = grey, midpoint = 0.5, name = expression(paste("E" [nl])), limits = c(0, 1), breaks =c(0, 0.5, 1));
p <- p + xlab("SNR") + ylab("Number of Cells") + labs(fill = "Localisation Error");
p <- p + scale_x_discrete(labels=c("0.10", "0.25", "0.50", "0.75", "1.00", "2.00", "3.00"));
p <- p + theme_minimal();
p + theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel);

#Cell centroid error
p <- ggplot(cellCountErrors, aes(as.character(snr), N_CELLS));
p <- p + geom_raster(aes(fill = MEAN_CELL_CENTROID_ERROR),interpolate=FALSE);
p <- p + scale_fill_gradient2(low=green, high=red, mid = yellow, midpoint = 6, name = expression(paste("E" [cl])), limits = c(1, 11), breaks =c(0, 5, 10));
p <- p + xlab("SNR") + ylab("Number of Cells") + labs(fill = "Localisation Error");
p <- p + scale_x_discrete(labels=c("0.10", "0.25", "0.50", "0.75", "1.00", "2.00", "3.00"));
p <- p + theme_minimal();
p + theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel);
