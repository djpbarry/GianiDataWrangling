allData$Inner <- allData$Distance_To_Centre < dist_thresh;

allData$VolRatio <- allData$Cell_Volume_microns / allData$Nucleus_Volume_microns;

allData$YapRatio <- allData$Nucleus_YAP1_Mean_Intensity / allData$Cytoplasm_YAP1_Mean_Intensity;

embryos <- unique(allData$Embryo);

cellData <- data.frame();

for(e in embryos){
  thisData <- allData[allData$Embryo == e,];
  allData[allData$Embryo == e, "n_Inner"] <- sum(thisData$Inner, na.rm = TRUE);
  allData[allData$Embryo == e, "n_Outer"] <- sum(!thisData$Inner, na.rm = TRUE);
  cellData <- rbind(cellData, data.frame(Inner = TRUE, nCells = sum(thisData$Inner, na.rm = TRUE), Treatment = thisData[1,treatmentHeading]));
  cellData <- rbind(cellData, data.frame(Inner = FALSE, nCells = sum(!thisData$Inner, na.rm = TRUE), Treatment = thisData[1,treatmentHeading]));
}

controlData <- subset(cellData, cellData$Treatment == CONTROL_VALUE);
treatedData <- subset(cellData, cellData$Treatment == TREATED_VALUE);

pdf(paste("plots", "inner_outer_cell_counts.pdf", sep=.Platform$file.sep), pointsize=20);

boxplot(nCells ~ Inner, data = controlData, outline = FALSE, main = "", ylab = "Number of cells", xlab = "", names = c("Outer", "Inner"), ylim=c(0,30));
beeswarm(nCells ~ Inner, data = controlData, pch = 16, cex = 0.7, col = c(2:4), add = TRUE);

dev.off();

controlData <- subset(allData, cellData$Treatment == CONTROL_VALUE);
treatedData <- subset(allData, cellData$Treatment == TREATED_VALUE);

pdf(paste("plots", "inner_outer_nuclear_volume.pdf", sep=.Platform$file.sep), pointsize=20);

boxplot(Nucleus_Volume_microns ~ Inner, data = controlData, outline = FALSE, main = "", ylab = "Nuclear Volume (microns^3)", xlab = "", names = c("Outer", "Inner"), ylim=c(0,100000));
beeswarm(Nucleus_Volume_microns ~ Inner, data = controlData, pch = 16, cex = 0.7, col = c(2:4), add = TRUE);

wilcox.test(controlData$Nucleus_Volume_microns ~controlData$Inner)

dev.off();

pdf(paste("plots", "inner_outer_cell_volume.pdf", sep=.Platform$file.sep), pointsize=20);

boxplot(Cell_Volume_microns ~ Inner, data = controlData, outline = FALSE, main = "", ylab = "Cell Volume (microns^3)", xlab = "", names = c("Outer", "Inner"), ylim=c(0,450000));
beeswarm(Cell_Volume_microns ~ Inner, data = controlData, pch = 16, cex = 0.7, col = c(2:4), add = TRUE);

wilcox.test(controlData$Cell_Volume_microns ~controlData$Inner)

dev.off();

pdf(paste("plots", "inner_outer_volume_ratio.pdf", sep=.Platform$file.sep), pointsize=20);

boxplot(VolRatio ~ Inner, data = controlData, outline = FALSE, main = "", ylab = "Nuclear/Cell Volume Ratio", xlab = "", names = c("Outer", "Inner"), ylim=c(0,15));
beeswarm(VolRatio ~ Inner, data = controlData, pch = 16, cex = 0.7, col = c(2:4), add = TRUE);

wilcox.test(controlData$VolRatio ~controlData$Inner)

dev.off();

pdf(paste("plots", "inner_outer_yap_ratio.pdf", sep=.Platform$file.sep), pointsize=20);

boxplot(YapRatio ~ Inner, data = controlData, outline = FALSE, main = "", ylab = "Nuclear/Cytoplasmic Yap1 Expression", xlab = "", names = c("Outer", "Inner"), ylim=c(0,8.0));
beeswarm(YapRatio ~ Inner, data = controlData, pch = 16, cex = 0.7, col = c(2:4), add = TRUE);

wilcox.test(controlData$YapRatio ~controlData$Inner)

dev.off();

#allData %>%
#  group_by(Embryo) %>%
#    summarise(
#      n = n(),
#      meanSize = mean(Cell_Volume_microns, na.rm = TRUE)
#    ) %>%
#      filter(n > 1)
