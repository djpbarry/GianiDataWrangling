source("published_scripts/Definitions.R");
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

#colours <- c("#FF0000", "#00DD00");

morphData <- allData;

morphData$Inner <- morphData$Distance_To_Centre < dist_thresh;

morphData$VolRatio <- morphData$Cell_Volume_microns / morphData$Nucleus_Volume_microns;

morphData$YapRatio <- morphData$Nucleus_YAP1_Mean_Intensity / morphData$Cytoplasm_YAP1_Mean_Intensity;

morphData$Description[morphData$Inner] <- "Inner";
morphData$Description[!morphData$Inner] <- "Outer";

embryos <- unique(morphData$Embryo);

cellData <- data.frame();

for(e in embryos){
  thisData <- morphData[morphData$Embryo == e,];
  morphData[morphData$Embryo == e, "n_Inner"] <- sum(thisData$Inner, na.rm = TRUE);
  morphData[morphData$Embryo == e, "n_Outer"] <- sum(!thisData$Inner, na.rm = TRUE);
  cellData <- rbind(cellData, data.frame(Inner = TRUE, nCells = sum(thisData$Inner, na.rm = TRUE), Treatment = thisData[1,treatmentHeading]));
  cellData <- rbind(cellData, data.frame(Inner = FALSE, nCells = sum(!thisData$Inner, na.rm = TRUE), Treatment = thisData[1,treatmentHeading]));
}

cellData$Description[cellData$Inner] <- "Inner";
cellData$Description[!cellData$Inner] <- "Outer";

controlData <- subset(cellData, cellData$Treatment == CONTROL_VALUE);
#treatedData <- subset(cellData, cellData$Treatment == TREATED_VALUE);

#pdf(paste("plots", "Fig3A.pdf", sep=.Platform$file.sep));

#boxplot(nCells ~ Inner, data = controlData, outline = FALSE, main = "", ylab = "Number of cells", xlab = "", names = c("Outer", "Inner"), ylim=c(0,35), lwd=1.0);
#beeswarm(nCells ~ Inner, data = controlData, pch = 16,cex=0.75, col = c(2:4), add = TRUE);

ggplot(controlData, aes(x=Description, y=nCells, color=Description)) + 
  geom_jitter(width=0.15,alpha=0.5,size=3) +
  geom_boxplot(width=0.6,color=black,coef=20,alpha=0.0) +
  scale_color_manual(values=c(blue,orange)) +
  #ylim(c(0,100)) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel, plot.title = axislabel) +
  guides(fill = "none", color="none") + xlab("") + ylab("Number of cells");

#dev.off();

controlData <- subset(morphData, morphData$Treatment == CONTROL_VALUE);
#treatedData <- subset(morphData, morphData$Treatment == TREATED_VALUE);

#pdf(paste("plots", "Fig3B.pdf", sep=.Platform$file.sep), width = 6, height = 6);

#boxplot(Nucleus_Volume_microns ~ Inner, data = controlData, outline = FALSE, main = "", ylab = expression("Nuclear Volume (" ~ mu ~ m^{3} ~ ")"), xlab = "", names = c("Outer", "Inner"), ylim=c(0,6000));
#beeswarm(Nucleus_Volume_microns ~ Inner, data = controlData, pch = 16, cex = 0.75, col = c(2:4), add = TRUE);
wilcox.test(controlData$Nucleus_Volume_microns ~controlData$Inner)
ks.test(controlData[ controlData$Inner,]$Nucleus_Volume_microns,controlData[!controlData$Inner,]$Nucleus_Volume_microns);

ggplot(controlData, aes(x=Description, y=Nucleus_Volume_microns, color=Description)) + 
  geom_jitter(width=0.15,alpha=0.5,size=3) +
  geom_boxplot(width=0.6,color=black,coef=20,alpha=0.0) +
  ylim(0,7000) +
  geom_segment(aes(x = 1, xend = 2, y=6000, yend=6000), color=black, size=0.75) +
  geom_text(label='n.s.',x=1.5, y=6400, color=black,size=6) +
  scale_color_manual(values=c(blue,orange)) +
  #ylim(c(0,100)) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel, plot.title = axislabel) +
  guides(fill = "none", color="none") + xlab("") + ylab(expression("Nuclear Volume (" ~ mu ~ m^{3} ~ ")"));

#dev.off();

#pdf(paste("plots", "inner_outer_cell_volume.pdf", sep=.Platform$file.sep), pointsize=20);

#boxplot(Cell_Volume_microns ~ Inner, data = controlData, outline = FALSE, main = "", ylab = expression("Cell Volume (" ~ mu ~ m^{3} ~ ")"), xlab = "", names = c("Outer", "Inner"), ylim = c(0,30000));
#beeswarm(Cell_Volume_microns ~ Inner, data = controlData, pch = 16, cex = 0.75, col = c(2:4), add = TRUE);

wilcox.test(controlData$Cell_Volume_microns ~controlData$Inner)

ggplot(controlData, aes(x=Description, y=Cell_Volume_microns, color=Description)) + 
  geom_jitter(width=0.15,alpha=0.5,size=3) +
  geom_boxplot(width=0.6,color=black,coef=20,alpha=0.0) +
  ylim(0,32000) +
  geom_segment(aes(x = 1, xend = 2, y=29000, yend=29000), color=black, size=0.75) +
  geom_text(label='n.s.',x=1.5, y=31000, color=black,size=6) +
  scale_color_manual(values=c(blue,orange)) +
  #ylim(c(0,100)) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel, plot.title = axislabel) +
  guides(fill = "none", color="none") + xlab("") + ylab(expression("Cell Volume (" ~ mu ~ m^{3} ~ ")"));



#dev.off();

#pdf(paste("plots", "inner_outer_volume_ratio.pdf", sep=.Platform$file.sep), pointsize=20);

#boxplot(VolRatio ~ Inner, data = controlData, outline = FALSE, main = "", ylab = "Cell/Nuclear Volume Ratio", xlab = "", names = c("Outer", "Inner"), ylim=c(1,200), log="y");
#beeswarm(VolRatio ~ Inner, data = controlData, pch = 16, cex = 0.75, col = c(2:4), add = TRUE);

wilcox.test(controlData$VolRatio ~controlData$Inner)


ggplot(controlData, aes(x=Description, y=VolRatio, color=Description)) + 
  geom_jitter(width=0.15,alpha=0.5,size=3) +
  geom_boxplot(width=0.6,color=black,coef=20,alpha=0.0) +
  geom_segment(aes(x = 1, xend = 2, y=200, yend=200), color=black, size=0.75) +
  geom_text(label='p < 0.001',x=1.5, y=2.5, color=black,size=6) +
  scale_color_manual(values=c(blue,orange)) +
  scale_y_log10("Cell/Nuclear Volume Ratio", limits=c(1,500)) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel, plot.title = axislabel) +
  guides(fill = "none", color="none") + xlab("") + ylab(expression("Cell/Nuclear Volume Ratio"));


#dev.off();

#pdf(paste("plots", "inner_outer_yap_ratio.pdf", sep=.Platform$file.sep), pointsize=20);

#boxplot(YapRatio ~ Inner, data = controlData, outline = FALSE, main = "", ylab = "Nuclear/Cytoplasmic Yap1 Expression", xlab = "", names = c("Outer", "Inner"), ylim=c(0.3,12), log="y");
#beeswarm(YapRatio ~ Inner, data = controlData, pch = 16, cex = 0.75, col = c(2:4), add = TRUE);

ggplot(controlData, aes(x=Description, y=YapRatio, color=Description)) + 
  geom_jitter(width=0.15,alpha=0.5,size=3) +
  geom_boxplot(width=0.6,color=black,coef=20,alpha=0.0) +
  geom_segment(aes(x = 1, xend = 2, y=12.5, yend=12.5), color=black, size=0.75) +
  geom_text(label='p < 0.0001',x=1.5, y=1.17, color=black,size=6) +
  scale_color_manual(values=c(blue,orange)) +
  scale_y_log10("Nuclear/Cytoplasmic YAP1 Expression", limits=c(0.3,15)) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel, plot.title = axislabel) +
  guides(fill = "none", color="none") + xlab("") + ylab(expression("Nuclear/Cytoplasmic YAP1 Expression"));

wilcox.test(controlData$YapRatio ~controlData$Inner)

#boxplot(YapRatio ~ Inner, data = treatedData, outline = FALSE, main = "", ylab = "Nuclear/Cytoplasmic Yap1 Expression", xlab = "", names = c("Outer", "Inner"), ylim=c(0.3,12), log="y");
#beeswarm(YapRatio ~ Inner, data = treatedData, pch = 16, cex = 0.75, col = c(2:4), add = TRUE);

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
