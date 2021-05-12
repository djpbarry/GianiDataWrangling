#library(Rmisc);
library(ggplot2);
#library(RColorBrewer);
#library(hrbrthemes);
#library(dplyr);

CONF_LEVEL <- 0.95;
gtDir <- "Z:/working/barryd/GIANI_Paper/Sim_Data/sim_images_ground_truth";
gianiDir <- "Z:/working/barryd/GIANI_Paper/Sim_Data/sim_images_output";
INDEX = "Index";
NUCLEUS <- "Nucleus";
CYTOPLASM <- "Cytoplasm";
CELL <- "Cell";
VOLUME_VOX <- "Volume_Vox";
VOLUME_MICRONS <- "Volume_microns";
SURFACE_AREA_VOXELS <- "Surface_Area_Voxels";
SURFACE_AREA_MICRONS <- "Surface_Area_microns";
GT_CELL_VOLUME_MIC <- "Cell_Volume_Microns_Cubed";
MEAN <- "Mean";
GIANI_FILE <- ".csv";
GROUND_TRUTH_FILE <- ".csv";
SNR <- "snr";
RUN <- "run_";
N_CELLS <- "ncells";
GT_CENTROID_X <- "Nucleus_Centroid_X";
GT_CENTROID_Y <- "Nucleus_Centroid_Y";
GT_CENTROID_Z <- "Nucleus_Centroid_Z";
GTC_CENTROID_X <- "Cell_Centroid_X";
GTC_CENTROID_Y <- "Cell_Centroid_Y";
GTC_CENTROID_Z <- "Cell_Centroid_Z";
GD_CENTROID_X <- "Centroid_X";
GD_CENTROID_Y <- "Centroid_Y";
GD_CENTROID_Z <- "Centroid_Z";
IM_CELL_VOL <- "Cell.Volume";
IM_PROP_VOL <- "Measured_Proportion_Of_Volume";
M_NUC_CENTROID_X <- "Measured_Nuc_Centroid_X";
M_NUC_CENTROID_Y <- "Measured_Nuc_Centroid_Y";
M_NUC_CENTROID_Z <- "Measured_Nuc_Centroid_Z";
M_CELL_CENTROID_X <- "Measured_Cell_Centroid_X";
M_CELL_CENTROID_Y <- "Measured_Cell_Centroid_Y";
M_CELL_CENTROID_Z <- "Measured_Cell_Centroid_Z";
GROUND_TRUTH_FOUND <- "Ground_Truth_Found";
NUC_CENTROID_ERROR <- "Nuc_Centroid_Error";
CELL_CENTROID_ERROR <- "Cell_Centroid_Error";
CELL_COUNT_ERROR <- "Cell_Count_Error";
PROP_CELL_COUNT_ERROR <- "Proportional_Cell_Count_Error";
LABEL <- "Label";
MEASURED_N <- "Measured_Cell_Count";
GROUND_TRUTH_N <- "Ground_Truth_Cell_Count";
VOL_ERROR <- "Volume_Error";
NORM_VOL_ERROR <- "Normalised_Volume_Error";
GT_PROP_VOL <- "Ground_Truth_Proportion_Of_Volume";
GD_PROP_VOL <- "Measured_Proportion_Of_Volume";
PROP_VOL_ERROR <- "Proportionate_Volume_Error";
NORM_DIST <- "Normalised_Distance_to_Centre";
MEAN_NUC_CENTROID_ERROR = "Mean_Nuclear_Centroid_Error";
MEAN_CELL_CENTROID_ERROR = "Mean_Cell_Centroid_Error";
YAP <- "YAP1";
GATA <- "GATA3";
DAPI <- "DAPI";
CONTROL_VALUE <- 1;
TREATED_VALUE <- 2;
treatmentHeading = "Treatment";
MEAN_INTENSITY <- "Mean_Intensity";
STD_INTENSITY <- "STD_Intensity";
MIN_INTENSITY <- "Min_Intensity";
MAX_INTENSITY <- "Max_Intensity";
INTEGRATED_DENSITY <- "Integrated_Density";

red <- "#DD0000";
green <- "#008800";
yellow <- "#FFFF00";
blue <- "#0000FF";
grey <- "#999999";
orange <- "#FF6600";
black <- "#000000";

dist_thresh <- 0.5;

cc <- scales::seq_gradient_pal(blue,orange,"Lab")(seq(0,1,length.out=4));

axislabel <- element_text(hjust=0.5, size=30, colour = "black");

euclidDistance <- function(v1, v2){
  if(length(v1) != length(v2)) return(NaN);
  sum <- 0.0;
  for(i in 1:length(v1)){
    sum <- sum + (v2[i] - v1[i])^2;
  }
  return(sqrt(sum))
}

buildDataFrame <- function(data_entries, nRows, nCols, label1, label2, headings){
  #browser();
  output <- data.frame(matrix(data=data_entries, nrow=nRows, ncol=nCols));

  colHeadings <- vector(mode="character",length=nCols);

  for(i in 1:nCols){
    if(is.null(label1)){
      colHeadings[i] <- headings[i];
    }else if(is.null(label2)){
      colHeadings[i] <- paste(label1, headings[i], sep="_");
    } else {
      colHeadings[i] <- paste(label1, label2, headings[i], sep="_");
    }
  }

  colnames(output) <- colHeadings;
  #browser();
  return(output);
}

saveHistogram <- function(outputDirectory, filename, controlData, treatedData, colour1, colour2, heading, xLabel, yLabel, xLimits, yLimits, binning){
  pdf(paste(outputDirectory, filename, sep=.Platform$file.sep));
  
  hist(controlData, col=colour1, border="black", main=heading, ylab=yLabel,xlab=xLabel, xlim=xLimits, ylim=yLimits, breaks=binning);
  hist(treatedData, col=colour2, border="black", breaks=binning, add=TRUE);
  legend("topright", c(CONTROL, TREATED), fill=c("red", "green"));
  
  dev.off();
}

saveConfidenceInterval <- function(data, headings, outputDir, filename){
  n <- ncol(data);
  
  browser();
  
  interval <- c(CI(data, CONF_LEVEL));
  
  ciResults <- data.frame(matrix(interval, nrow=1,ncol=3 * n));
  
  ciHeadings <- c(paste(headings[1], "upper", CONF_LEVEL * 100, "% limit"),
                  paste(headings[1], "mean"),
                  paste(headings[1], "lower", CONF_LEVEL * 100, "% limit"));
  
  for(i in 2:n){
    ciHeadings <- c(paste(headings[i], "upper", CONF_LEVEL * 100, "% limit"),
                             paste(headings[i], "mean"),
                             paste(headings[i], "lower", CONF_LEVEL * 100, "% limit"));
  }
  
  colnames(ciResults) <- ciHeadings;
  
  write.csv(ciResults, paste(outputDir, filename, sep="/"));
}

saveConfidenceIntervals <- function(controlData, treatedData, outputDir, filename){
  interval <- c(CI(controlData, CONF_LEVEL), CI(treatedData, CONF_LEVEL));
  
  ciResults <- data.frame(matrix(interval, nrow=1,ncol=6));
  colnames(ciResults) <- c(paste(CONTROL, "upper", CONF_LEVEL * 100, "% limit"),
                                  paste(CONTROL, "mean"),
                                  paste(CONTROL, "lower", CONF_LEVEL * 100, "% limit"),
                                  paste(TREATED, "upper", CONF_LEVEL * 100, "% limit"),
                                  paste(TREATED, "mean"),
                                  paste(TREATED, "lower", CONF_LEVEL * 100, "% limit"));
  
  write.csv(ciResults, paste(outputDir, filename, sep="/"));
}

savePlot <- function(outputDir, filename, x, controlData, treatedData, controlErr, treatedErr, title, xLabel, yLabel, yLimits){
  pdf(paste(outputDir, filename, sep=.Platform$file.sep));
  
  plot(x, controlData, xlab=xLabel, main=title, ylab=yLabel, ylim=yLimits, col="red", pch=15);
  arrows(x, controlData-controlErr, x, controlData+controlErr, length=0.05, angle=90, code=3, col="red");
  
  points(x, treatedData, col="green", pch=15);
  arrows(x, treatedData-treatedErr, x, treatedData+treatedErr, length=0.05, angle=90, code=3, col="green");
  
  legend("topright", c(CONTROL, TREATED), fill=c("red", "green"));
  
  dev.off();
}

getMeansErrors <- function(params, conditions, inputData){
  
  length <- length(inputData[[j]][[t]]);
  means <- matrix(nrow=length,ncol=params*conditions);
  err <- matrix(nrow=length,ncol=params*conditions);
  
  for(j in 1:params){
    for(t in 1:conditions){
      for(i in 1:length){
        if(length(inputData[[j]][[t]][[i]]) > 1){
          means[i, t + (j-1)*2] <- mean(unlist(inputData[[j]][[t]][[i]]));
          interval <- CI(unlist(inputData[[j]][[t]][[i]]), CONF_LEVEL);
          err[i,t+(j-1)*2] <- interval[1] - interval[2];
        } else {
          means[i,t+(j-1)*2] <- unlist(inputData[[j]][[t]][[i]]);
          err[i,t+(j-1)*2] <- 0.0;
        }
      }
    }
  }
  return(list(means, err));
}

saveMeanErrors <- function(x, means, err, labels, outputDir, filename){
  meanErrDataFrame <- data.frame(cbind(x, means, err));
  newColLabels <- c(labels[1]);
  for(i in 2:length(labels)){
    newColLabels <- cbind(newColLabels,
                          paste(CONTROL, MEAN, labels[i]),
                          paste(TREATED, MEAN, labels[i]));
  }
  for(i in 2:length(labels)){
    newColLabels <- cbind(newColLabels,
                          paste(CONTROL, CONF_LEVEL*100, "% CI of ", labels[i]),
                          paste(TREATED, CONF_LEVEL*100, "% CI of ", labels[i]));
  }
  
  colnames(meanErrDataFrame) <- newColLabels;
  
  write.csv(meanErrDataFrame, paste(outputDir, filename, sep="/"));
}
