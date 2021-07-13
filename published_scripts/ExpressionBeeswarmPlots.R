source("published_scripts/definitions.R");

par(cex=1.5,
    cex.axis=1.5,
    cex.lab = 1.5,mai=c(0.75,1.5,0.5,0.5));

expData <- allData;

expData$Inner[expData$Distance_To_Centre < dist_thresh] <- "Inner";
expData$Inner[expData$Distance_To_Centre >= dist_thresh] <- "Outer";

expData$NormGata <- expData$Nucleus_GATA3_Mean_Intensity / expData$Nucleus_DAPI_Mean_Intensity;
expData$NormYap <- expData$Nucleus_YAP1_Mean_Intensity/expData$Cytoplasm_YAP1_Mean_Intensity;

expData$Description[expData$Treatment==CONTROL_VALUE] <- paste(expData$Inner[expData$Treatment==CONTROL_VALUE], "Control");
expData$Description[expData$Treatment==TREATED_VALUE] <- paste(expData$Inner[expData$Treatment==TREATED_VALUE], "Treated");

innerControl <- expData[ expData$Description == "Inner Control",];
innerControl <- innerControl[order(-innerControl$NormYap),];

outerControl <- expData[ expData$Description == "Outer Control",];
outerControl <- outerControl[order(outerControl$NormYap),];

#pdf(paste("plots", "control_treated_yap1.pdf", sep=.Platform$file.sep), pointsize=13);

gataData <- expData[complete.cases(expData),];

#boxplot(NormGata ~ Description, data = gataData, outline = FALSE, main = "", ylim = c(0.02,1.5), ylab = "Nuclear GATA3/DAPI Expression", xlab = "", names = c("Outer_Control", "Outer_Treated", "Inner_Control", "Inner_Treated"), log="y");
#beeswarm(NormGata ~ Description, data = gataData, pch = 16, cex = 0.75, col = c(2:5), add = TRUE);

ggplot(gataData, aes(x=Description, y=NormGata, color=Description)) + 
  geom_jitter(width=0.15,alpha=0.5,size=3) +
  geom_boxplot(width=0.6,color=black,coef=20,alpha=0.0) +
  geom_segment(aes(x = 1, xend = 3, y=2.1, yend=2.1), color=black, size=0.75) +
  geom_text(label='p < 0.001',x=2, y=0.42, color=black,size=10) +
  geom_segment(aes(x = 2, xend = 4, y=3.4, yend=3.4), color=black, size=0.75) +
  geom_text(label='n.s.',x=3, y=0.61, color=black,size=10) +
  geom_segment(aes(x = 1, xend = 2, y=1.4, yend=1.4), color=black, size=0.75) +
  geom_text(label='n.s.',x=1.5, y=0.24, color=black,size=10) +
  geom_segment(aes(x = 3, xend = 4, y=1.4, yend=1.4), color=black, size=0.75) +
  geom_text(label='p < 0.01',x=3.5, y=0.24, color=black,size=10) +
  scale_color_manual(values=cc2) +
  #ylim(c(0,100)) +
  theme_linedraw() +
  scale_y_log10("Nuclear GATA3/DAPI\n Expression", limits=c(0.02,4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel, plot.title = axislabel) +
  guides(fill = "none", color="none") + xlab("") + 
  scale_x_discrete(breaks=unique(gataData$Description), labels = c("Inner\nControl", "Inner\nTreated", "Outer\nControl", "Outer\nTreated"));




wilcox.test(subset(gataData, gataData$Treatment == 1)$NormGata ~ subset(gataData, gataData$Treatment == 1)$Inner)
wilcox.test(subset(gataData, gataData$Treatment == 2)$NormGata ~ subset(gataData, gataData$Treatment == 2)$Inner)

wilcox.test(subset(gataData, gataData$Inner == "Outer")$NormGata ~ subset(gataData, gataData$Inner == "Outer")$Treatment)
wilcox.test(subset(gataData, gataData$Inner == "Inner")$NormGata ~ subset(gataData, gataData$Inner == "Inner")$Treatment)

?subset

#dev.off();

#pdf(paste("plots", "control_treated_gata3.pdf", sep=.Platform$file.sep), pointsize=13);

#boxplot(NormYap ~ Description, data = expData, log = "y", outline = FALSE, ylim = c(0.3,15.0), main = "", ylab = "Nuclear/Cytoplasmic YAP1 Expression", xlab = "", names = c("Outer_Control", "Outer_Treated", "Inner_Control", "Inner_Treated"));
#beeswarm(NormYap ~ Description, data = expData, pch = 16, cex = 0.7, col = c(2:4), add = TRUE);


ggplot(expData, aes(x=Description, y=NormYap, color=Description)) + 
  geom_jitter(width=0.15,alpha=0.5,size=3) +
  geom_boxplot(width=0.6,color=black,coef=20,alpha=0.0) +
  geom_segment(aes(x = 1, xend = 3, y=16, yend=16), color=black, size=0.75) +
  geom_text(label='p < 0.0001',x=2, y=1.28, color=black,size=10) +
  geom_segment(aes(x = 2, xend = 4, y=23, yend=23), color=black, size=0.75) +
  geom_text(label='p < 0.0001',x=3, y=1.43, color=black,size=10) +
  geom_segment(aes(x = 1, xend = 2, y=12.5, yend=12.5), color=black, size=0.75) +
  geom_text(label='n.s.',x=1.5, y=1.17, color=black,size=10) +
  geom_segment(aes(x = 3, xend = 4, y=12.5, yend=12.5), color=black, size=0.75) +
  geom_text(label='p < 0.0001',x=3.5, y=1.17, color=black,size=10) +
  scale_color_manual(values=cc2) +
  #ylim(c(0,100)) +
  theme_linedraw() +
  scale_y_log10("Nuclear/Cytoplasmic\nYAP1 Expression", limits=c(0.2,25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel, plot.title = axislabel) +
  guides(fill = "none", color="none") + xlab("") +
  scale_x_discrete(breaks=unique(gataData$Description), labels = c("Inner\nControl", "Inner\nTreated", "Outer\nControl", "Outer\nTreated"));




wilcox.test(subset(expData, expData$Inner == "Outer")$NormYap ~ subset(expData, expData$Inner == "Outer")$Treatment)
wilcox.test(subset(expData, expData$Inner == "Inner")$NormYap ~ subset(expData, expData$Inner == "Inner")$Treatment)

wilcox.test(subset(expData, expData$Treatment == 1)$NormYap ~ subset(expData, expData$Treatment == 1)$Inner)
wilcox.test(subset(expData, expData$Treatment == 2)$NormYap ~ subset(expData, expData$Treatment == 2)$Inner)


#dev.off();

