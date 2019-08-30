source("exp_scripts/Definitions.R");

nuc <- 1;
cyto <- 3;
cell <- 5;

files <- list.files(path = directory, pattern = ".csv", recursive = TRUE, full.names = TRUE);

allData <- list(data.frame(x=NULL), data.frame(x=NULL), data.frame(x=NULL), data.frame(x=NULL), data.frame(x=NULL), data.frame(x=NULL));

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

channelLabels <- c(GATA, YAP);

embryoIndex <- 1;

for (f in files){
  if(grepl("18.06.10", f) || (grepl("19.06.10", f) && !(grepl("_S13_", f) || grepl("_S15_", f) || grepl("_S6_", f)))){
    channel <- abs(channel - 1);
  }
  thisData <- read.csv(f,  encoding="UTF-8", stringsAsFactors=FALSE);
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
    if(thisData$Channel[i] == 1){
      channel = 0;
    } else{
      channel = 1;
    }
    if(grepl("18.06.10", f) || (grepl("19.06.10", f) && !(grepl("_S13_", f) || grepl("_S15_", f) || grepl("_S6_", f)))){
      channel <- abs(channel - 1);
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
    intensData <- buildDataFrame(c(thisData$Mean.Pixel.Value[i],
                                   thisData$Pixel.Standard.Deviation[i],
                                   thisData$Min.Pixel.Value[i],
                                   thisData$Max.Pixel.Value[i],
                                   thisData$Integrated.Density[i]),
                                 1,
                                 5,
                                 regionLabel,
                                 channelLabels[channel + 1],
                                 c(MEAN_INTENSITY,
                                   STD_INTENSITY,
                                   MIN_INTENSITY,
                                   MAX_INTENSITY,
                                   INTEGRATED_DENSITY));
    morphData <- NULL;
    if(channel == 0){
      dataEntries <- c(thisData$Volume..Voxels.[i], thisData$Volume...U.00B5.m.3.[i], thisData$Surface.Area..Voxels.[i], thisData$Surface.Area...U.00B5.m.2.[i]);
      morphData <- buildDataFrame(dataEntries,
                                  1,
                                  4,
                                  regionLabel,
                                  NULL,
                                  c(VOLUME_VOX,
                                    VOLUME_MICRONS,
                                    SURFACE_AREA_VOXELS,
                                    SURFACE_AREA_MICRONS));
    }
    if(channel == 0){
      if(identical(regionLabel, NUCLEUS)){
        row <- cbind(data.frame(Label=f,
                                Embryo=embryoIndex,
                                Index=thisData$Index[i],
                                Distance_To_Centre=thisData$Normalised.Distance.to.Centre[i]),
                     intensData,
                     morphData,
                     treated);
      } else {
        row <- cbind(intensData, morphData);
      }
    } else {
      row <- intensData;
    }
    frameIndex <- label + channel;
    allData[[frameIndex]] <- rbind(allData[[frameIndex]], row);
  }
  embryoIndex <- embryoIndex + 1;
}

allData <- do.call(cbind, allData);
