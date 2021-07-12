source("exp_scripts/Definitions.R");
#c(bottom, left, top, right)

par(cex=1.5,
  cex.axis=1.5,
    cex.lab = 1.5,mai=c(0.75,1.5,0.5,0.5));



#xAxislabel <- element_text(vjust=-1, size=35, colour = "black");
#yAxislabel <- element_text(size=35, colour = "black");

#customTheme <- theme(axis.text.y = yAxislabel,
#                     axis.text.x = xAxislabel,
#                     axis.title.x=xAxislabel,
#                     axis.title.y = yAxislabel,
#                     axis.ticks = element_line(colour = "black", size = 1),
#                     axis.ticks.length=unit(.4, "cm"),
#                     legend.position = "none",
#                     panel.border = element_rect(colour = "black", fill=NA, size=1));

colours <- c("#FF0000", "#00DD00");

morphData <- allData;

morphData$Inner <- morphData$Distance_To_Centre < dist_thresh;

morphData$VolRatio <- morphData$Cell_Volume_microns / morphData$Nucleus_Volume_microns;

morphData$YapRatio <- morphData$Nucleus_YAP1_Mean_Intensity / morphData$Cytoplasm_YAP1_Mean_Intensity;

embryos <- unique(morphData$Embryo);

cellData <- data.frame();

for(e in embryos){
  thisData <- morphData[morphData$Embryo == e,];
  morphData[morphData$Embryo == e, "n_Inner"] <- sum(thisData$Inner, na.rm = TRUE);
  morphData[morphData$Embryo == e, "n_Outer"] <- sum(!thisData$Inner, na.rm = TRUE);
  cellData <- rbind(cellData, data.frame(Inner = TRUE, nCells = sum(thisData$Inner, na.rm = TRUE), Treatment = thisData[1,treatmentHeading]));
  cellData <- rbind(cellData, data.frame(Inner = FALSE, nCells = sum(!thisData$Inner, na.rm = TRUE), Treatment = thisData[1,treatmentHeading]));
}

controlData <- subset(cellData, cellData$Treatment == CONTROL_VALUE);
treatedData <- subset(cellData, cellData$Treatment == TREATED_VALUE);

#pdf(paste("plots", "Fig3A.pdf", sep=.Platform$file.sep));

boxplot(nCells ~ Inner, data = controlData, outline = FALSE, main = "", ylab = "Number of cells", xlab = "", names = c("Outer", "Inner"), ylim=c(0,35), lwd=1.0);
beeswarm(nCells ~ Inner, data = controlData, pch = 16,cex=0.75, col = c(2:4), add = TRUE);

# p <- ggplot(cellData, aes(Inner, nCells, fill=Inner, color=Inner));
# p <- p + stat_boxplot(geom = "errorbar", color = "black", width = 0.5, lwd = 1);
# p <- p + geom_boxplot(outlier.colour = "white", notch = TRUE, color="black", lwd = 1, fill="white");
# #p <- p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1);
# p <- p + geom_jitter(shape=21, position=position_jitter(0.05), size=3);
# p <- p + scale_fill_manual(values=colours);
# p <- p + scale_color_manual(values=colours);
# p <- p + ylab("Number of Cells") + xlab("");
# p <- p + theme_classic();
# p + customTheme;

#dev.off();

controlData <- subset(morphData, morphData$Treatment == CONTROL_VALUE);
treatedData <- subset(morphData, morphData$Treatment == TREATED_VALUE);

#pdf(paste("plots", "Fig3B.pdf", sep=.Platform$file.sep), width = 6, height = 6);

boxplot(Nucleus_Volume_microns ~ Inner, data = controlData, outline = FALSE, main = "", ylab = expression("Nuclear Volume (" ~ mu ~ m^{3} ~ ")"), xlab = "", names = c("Outer", "Inner"), ylim=c(0,6000));
beeswarm(Nucleus_Volume_microns ~ Inner, data = controlData, pch = 16, cex = 0.75, col = c(2:4), add = TRUE);
wilcox.test(controlData$Nucleus_Volume_microns ~controlData$Inner)

# p <- ggplot(morphData, aes(Inner, Nucleus_Volume_microns, fill=Inner, color=Inner));
# p <- p + stat_boxplot(geom = "errorbar", color = "black", width = 0.5, lwd = 1);
# p <- p + geom_boxplot(outlier.colour = "white", notch = TRUE, color="black", lwd = 1, fill="white");
# p <- p + geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position = position_jitter(width = 0.05, height=0.05));
# p <- p + scale_fill_manual(values=colours);
# p <- p + scale_color_manual(values=colours);
# p <- p + ylab("Number of Cells") + xlab("");
# p <- p + theme_classic();
# p + customTheme;

#dev.off();

#pdf(paste("plots", "inner_outer_cell_volume.pdf", sep=.Platform$file.sep), pointsize=20);

boxplot(Cell_Volume_microns ~ Inner, data = controlData, outline = FALSE, main = "", ylab = expression("Cell Volume (" ~ mu ~ m^{3} ~ ")"), xlab = "", names = c("Outer", "Inner"), ylim = c(0,30000));
beeswarm(Cell_Volume_microns ~ Inner, data = controlData, pch = 16, cex = 0.75, col = c(2:4), add = TRUE);

wilcox.test(controlData$Cell_Volume_microns ~controlData$Inner)

#dev.off();

#pdf(paste("plots", "inner_outer_volume_ratio.pdf", sep=.Platform$file.sep), pointsize=20);

boxplot(VolRatio ~ Inner, data = controlData, outline = FALSE, main = "", ylab = "Cell/Nuclear Volume Ratio", xlab = "", names = c("Outer", "Inner"), ylim=c(1,200), log="y");
beeswarm(VolRatio ~ Inner, data = controlData, pch = 16, cex = 0.75, col = c(2:4), add = TRUE);

wilcox.test(controlData$VolRatio ~controlData$Inner)

#dev.off();

#pdf(paste("plots", "inner_outer_yap_ratio.pdf", sep=.Platform$file.sep), pointsize=20);

boxplot(YapRatio ~ Inner, data = controlData, outline = FALSE, main = "", ylab = "Nuclear/Cytoplasmic Yap1 Expression", xlab = "", names = c("Outer", "Inner"), ylim=c(0.3,12), log="y");
beeswarm(YapRatio ~ Inner, data = controlData, pch = 16, cex = 0.75, col = c(2:4), add = TRUE);

wilcox.test(controlData$YapRatio ~controlData$Inner)

boxplot(YapRatio ~ Inner, data = treatedData, outline = FALSE, main = "", ylab = "Nuclear/Cytoplasmic Yap1 Expression", xlab = "", names = c("Outer", "Inner"), ylim=c(0.3,12), log="y");
beeswarm(YapRatio ~ Inner, data = treatedData, pch = 16, cex = 0.75, col = c(2:4), add = TRUE);

wilcox.test(treatedData$YapRatio ~treatedData$Inner)

#dev.off();

#morphData %>%
#  group_by(Embryo) %>%
#    summarise(
#      n = n(),
#      meanSize = mean(Cell_Volume_microns, na.rm = TRUE)
#    ) %>%
#      filter(n > 1)

temp <- subset(morphData, select=c("Distance_To_Centre", "Inner"));
