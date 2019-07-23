source("scripts/Definitions.R");

embryoData <- allData[[embryoHeading]];
embryos <- unique(embryoData);
cellMorph <- list(data.frame(x=NULL), data.frame(x=NULL));

for(e in embryos){
  thisData <- subset(allData, embryoData==e);
  treatmentData <- thisData[[treatmentHeading]];
  if(treatmentData[1] > CONTROL_VALUE){
    label = TREATED;
  } else {
    label = CONTROL;
  }
  row <- morphData <- buildDataFrame(c(mean(thisData[[NUCLEUS_VOLUME_MIC]]),
                                       mean(thisData[[CELL_VOLUME_MIC]]),
                                       mean(thisData[[NUCLEUS_VOLUME_MIC]] / thisData[[CELL_VOLUME_MIC]])),
                                     1,
                                     3,
                                     label,
                                     NULL, c(paste("Mean", NUCLEUS_VOLUME_MIC, sep="_"),
                                             paste("Mean", CELL_VOLUME_MIC, sep="_"),
                                             paste("Mean", NUC_TO_CELL_VOLUME_RATIO, sep="_")));
  cellMorph[[treatmentData[1]]] <- rbind(cellMorph[[treatmentData[1]]], row);
}

columns <- c(colnames(cellMorph[[CONTROL_VALUE]]), colnames(cellMorph[[TREATED_VALUE]]));

pdf(paste("plots", "nuc_vol.pdf", sep=.Platform$file.sep));

hist(cellMorph[[CONTROL_VALUE]][[columns[1]]], col=rgb(0.8,0.0,0.0,0.5), border="black", main="Nuclear Volume", ylab="Number of Embryos",xlab=expression(paste("Mean Nuclear Volume (",mu,"m^3)")), xlim=c(1000,3500), ylim=c(0,10), breaks=seq(from=1000, to=3500, by=250));
hist(cellMorph[[TREATED_VALUE]][[columns[4]]], col=rgb(0.0,0.8,0.0,0.5), border="black", breaks=seq(from=1000, to=3500, by=250), add=TRUE);
legend("topright", c(CONTROL, TREATED), fill=c("red", "green"));

dev.off();

pdf(paste("plots", "cell_vol.pdf", sep=.Platform$file.sep));

hist(cellMorph[[CONTROL_VALUE]][[columns[2]]], col=rgb(0.8,0.0,0.0,0.5), border="black", main="Cell Volume", ylab="Number of Embryos",xlab=expression(paste("Mean Cell Volume (",mu,"m^3)")), xlim=c(4000,18000), ylim=c(0,8), breaks=seq(from=4000, to=18000, by=2000));
hist(cellMorph[[TREATED_VALUE]][[columns[5]]], col=rgb(0.0,0.8,0.0,0.5), border="black", breaks=seq(from=4000, to=18000, by=2000), add=TRUE);
legend("topright", c(CONTROL, TREATED), fill=c("red", "green"));

dev.off();

pdf(paste("plots", "nuc_to_cell_vol_ratio.pdf", sep=.Platform$file.sep));

hist(cellMorph[[CONTROL_VALUE]][[columns[3]]], col=rgb(0.8,0.0,0.0,0.5), border="black", main="Nuclear/Cell Volume Ratio", ylab="Number of Embryos",xlab="Mean Nuclear:Cell Volume Ratio", xlim=c(0.1,0.5), ylim=c(0,10), breaks=seq(from=0.1, to=0.5, by=0.05));
hist(cellMorph[[TREATED_VALUE]][[columns[6]]], col=rgb(0.0,0.8,0.0,0.5), border="black", breaks=seq(from=0.1, to=0.5, by=0.05), add=TRUE);
legend("topright", c(CONTROL, TREATED), fill=c("red", "green"));

dev.off();

#cellMorph <- merge(cellMorph[[CONTROL_VALUE]], cellMorph[[TREATED_VALUE]], by = 0, all = TRUE)[-1];

#write.csv(cellMorph, paste("outputs", "cell_morphology.csv", sep="/"));