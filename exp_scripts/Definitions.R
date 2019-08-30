library(Rmisc);

CONF_LEVEL <- 0.95;
directory <- "D:/OneDrive - The Francis Crick Institute/Working Data/Niakan/Claudia/GIANI_Paper";
treatmentHeading = "Treatment";
CONTROL_VALUE <- 1;
TREATED_VALUE <- 2;
embryoHeading = "Embryo";
NUCLEUS <- "Nucleus";
CYTOPLASM <- "Cytoplasm";
CELL <- "Cell";
YAP <- "YAP1";
GATA <- "GATA3";
MEAN_INTENSITY <- "Mean_Intensity";
STD_INTENSITY <- "STD_Intensity";
MIN_INTENSITY <- "Min_Intensity";
MAX_INTENSITY <- "Max_Intensity";
INTEGRATED_DENSITY <- "Integrated_Density";
VOLUME_VOX <- "Volume_Vox";
VOLUME_MICRONS <- "Volume_microns";
SURFACE_AREA_VOXELS <- "Surface_Area_Voxels";
SURFACE_AREA_MICRONS <- "Surface_Area_microns";
NUCLEUS_SURFACE_AREA_MIC <- paste(NUCLEUS, SURFACE_AREA_MICRONS, sep="_");
NUCLEUS_VOLUME_MIC <- paste(NUCLEUS, VOLUME_MICRONS, sep="_");
CELL_SURFACE_AREA_MIC <- paste(CELL, SURFACE_AREA_MICRONS, sep="_");
CELL_VOLUME_MIC <- paste(CELL, VOLUME_MICRONS, sep="_");
DISTANCE_TO_CENTRE <- "Distance_To_Centre";
NUC_TO_CELL_VOLUME_RATIO <- paste(NUCLEUS, "to", CELL, "Volume Ratio", sep="_");
VERSUS_DISTANCE <- paste("Versus", "Distance", sep="_");
CELL_VOLUME_MICRONS_VERSUS_DISTANCE <- paste(CELL_VOLUME_MIC, VERSUS_DISTANCE, sep="_");
GATA_VERSUS_DISTANCE <- paste(GATA, VERSUS_DISTANCE, sep="_");
YAP_VERSUS_DISTANCE <- paste(YAP, "Ratio", VERSUS_DISTANCE, sep="_");

TREATED <- "Treated";
CONTROL <- "Control";
MEAN <- "Mean";

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
  pdf(paste(outputDirectory, filename, sep=.Platform$file.sep), pointsize=20);
  
  hist(controlData, col=colour1, border="black", main="", ylab=yLabel,xlab=xLabel, xlim=xLimits, ylim=yLimits, breaks=binning);
  hist(treatedData, col=colour2, border="black", breaks=binning, add=TRUE);
  legend("topright", c(CONTROL, TREATED), fill=c("red", "green"));
  
  dev.off();
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
  pdf(paste(outputDir, filename, sep=.Platform$file.sep), pointsize=20);
  
  plot(x, controlData, xlab=xLabel, main="", ylab=yLabel, ylim=yLimits, col="red", pch=15);
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
