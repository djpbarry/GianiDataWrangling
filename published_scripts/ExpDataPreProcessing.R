source("published_scripts/Definitions.R");

nuc <- 1;
cyto <- 2;
cell <- 3;

directory <- "Z:/working/barryd//GIANI_Paper/experimental_data";

files <- list.files(path = directory, pattern = ".csv", recursive = TRUE, full.names = TRUE);

allData <- vector(mode = "list", length = 3);

treatedPatterns <- c("/18.06.09LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S4_Output",
                     "/18.06.09LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S7_Output",
                     "/18.06.09LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S8_Output",
                     "/18.07.22/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S1_Output",
                     "/18.07.22/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S6_Output",
                     "/18.07.22/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S9_Output",
                     "/18.07.22/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S13_Output",
                     "/18.07.22/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S15_Output",
                     "/18.07.22/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S16_Output",
                     "/18.06.10/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S10_Output",
                     "/18.06.09/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S4_Output",
                     "/18.06.09/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S7_Output",
                     "/18.06.09/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S8_Output",
                     "/18.06.25LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S7_Output",
                     "/18.06.25LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S11_Output",
                     "/19.06.10LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S17_Output",
                     "/19.06.10LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S18_Output",
                     "/19.06.10LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S19_Output",
                     "/19.06.10LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S20_Output",
                     "/19.06.10LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S13_Output",
                     "/19.06.10LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S15_Output",
                     "/19.07.05LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S7_Output",
                     "/19.07.05LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S9_Output");

controlPatterns <- c("/18.07.22/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S19_Output",
                     "/18.07.22/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S20_Output",
                     "/18.07.22/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S22_Output",
                     "/18.07.22/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S24_Output",
                     "/18.06.10/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S5_Output",
                     "/18.06.10/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S8_Output",
                     "/18.06.10/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S9_Output",
                     "/18.06.25LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S15_Output",
                     "/18.06.25LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S14_Output",
                     "/18.06.25LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S3_Output",
                     "/19.06.10LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S26_Output",
                     "/19.06.10LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S27_Output",
                     "/19.06.10LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S28_Output",
                     "/19.06.10LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S30_Output",
                     "/19.06.10LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S6_Output",
                     "/19.07.05LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S11_Output",
                     "/19.07.05LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S13_Output",
                     "/19.07.05LeicaSp5/GIANI v2.[[:digit:]]{3}_[[:print:]]+_S14_Output");

channelLabels <- c(DAPI, GATA, YAP);
channelIndices <- c(0, 1, 4);

embryoIndex <- 1;

for (f in files){
  thisData <- read.csv(f,  encoding="UTF-16", stringsAsFactors=FALSE);
  cols <- colnames(thisData);
  n <- nrow(thisData);
  for(tp in treatedPatterns){
    if(grepl(tp, f)){
      treated <- data.frame(Treatment=TREATED_VALUE);
    }
  }
  for(cp in controlPatterns){
    if(grepl(cp, f)){
      treated <- data.frame(Treatment=CONTROL_VALUE);
    }
  }
  for(i in 1:n){
    if(grepl("18.06.10", f) || (grepl("19.06.10", f) && !(grepl("_S13_", f) || grepl("_S15_", f) || grepl("_S6_", f)))){
      reverse <- TRUE;
    } else {
      reverse <- FALSE;
    }
    if(grepl(NUCLEUS, thisData$Label[i])){
      label <- nuc;
      regionLabel <- NUCLEUS;
    } else if(grepl(CYTOPLASM, thisData$Label[i])){
      label <- cyto;
      regionLabel <- CYTOPLASM;
    } else{
      label <- cell;
      regionLabel <- CELL;
    }
    intensData <- vector(mode = "list", length = 3);
    for(c in 1:3){
      channelName = channelLabels[c];
      if(c == 2 && reverse){
        channelName = channelLabels[3];
      } else if(c == 3 && reverse){
        channelName = channelLabels[2];
      }
      if(grepl("19.07.05", f) && c == 2){
        intensData[[c]] <- buildDataFrame(rep(NaN, 5),
                                          1,
                                          5,
                                          regionLabel,
                                          channelName,
                                          c(MEAN_INTENSITY,
                                            STD_INTENSITY,
                                            MIN_INTENSITY,
                                            MAX_INTENSITY,
                                            INTEGRATED_DENSITY));
      } else {
        intensData[[c]] <- buildDataFrame(c(thisData[[paste("Mean_Pixel_Value_C", channelIndices[c], sep="")]][i],
                                            thisData[[paste("Pixel_Standard_Deviation_C", channelIndices[c], sep="")]][i],
                                            thisData[[paste("Min_Pixel_Value_C", channelIndices[c], sep="")]][i],
                                            thisData[[paste("Max_Pixel_Value_C", channelIndices[c], sep="")]][i],
                                            thisData[[paste("Integrated_Density_C", channelIndices[c], sep="")]][i]),
                                     1,
                                     5,
                                     regionLabel,
                                     channelName,
                                     c(MEAN_INTENSITY,
                                       STD_INTENSITY,
                                       MIN_INTENSITY,
                                       MAX_INTENSITY,
                                       INTEGRATED_DENSITY));
      }
    }
    intensData <- do.call(cbind, intensData);
    dataEntries <- c(thisData$Volume..Voxels.[i], thisData$Volume..µm.3.[i], thisData$Surface.Area..Voxels.[i], thisData$Surface.Area..µm.2.[i]);
    morphData <- buildDataFrame(dataEntries,
                                1,
                                4,
                                regionLabel,
                                NULL,
                                c(VOLUME_VOX,
                                  VOLUME_MICRONS,
                                  SURFACE_AREA_VOXELS,
                                  SURFACE_AREA_MICRONS));
    if(identical(regionLabel, NUCLEUS)){
      row <- cbind(data.frame(Label=f,
                              Embryo=embryoIndex,
                              Index=thisData$Index[i],
                              Distance_To_Centre=thisData$Normalised_Distance_to_Centre[i]),
                   intensData,
                   morphData,
                   treated);
    } else {
      row <- cbind(data.frame(Embryo=embryoIndex, Index=thisData$Index[i]), intensData, morphData);
    }
    allData[[label]] <- rbind(allData[[label]], row);
  }
  #print(paste(embryoIndex, nrow(allData[[1]]), nrow(allData[[2]]), nrow(allData[[3]])));
  embryoIndex <- embryoIndex + 1;
}

#allData <- do.call(cbind, allData);

allData1 <- merge(allData[[1]], allData[[2]], by=c("Embryo", "Index"));
allData <- merge(allData1, allData[[3]], by=c("Embryo", "Index"));
