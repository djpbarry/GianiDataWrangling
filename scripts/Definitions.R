library(Rmisc);

confLevel <- 0.95;
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

  return(output);
}