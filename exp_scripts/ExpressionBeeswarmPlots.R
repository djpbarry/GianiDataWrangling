allData$NormGata <- allData$Nucleus_GATA3_Mean_Intensity / allData$Nucleus_DAPI_Mean_Intensity;
allData$NormYap <- allData$Nucleus_YAP1_Mean_Intensity / allData$Nucleus_DAPI_Mean_Intensity;

allData$Description <- paste(allData$Inner, allData$Treatment, sep = "_");

pdf(paste("plots", "control_treated_yap1.pdf", sep=.Platform$file.sep), pointsize=13);

boxplot(NormGata ~ Description, data = allData, log = "y", outline = FALSE, main = "", ylim = c(0.01,2.0), ylab = "Nuclear GATA3/DAPI Expression", xlab = "", names = c("Outer_Control", "Outer_Treated", "Inner_Control", "Inner_Treated"));
beeswarm(NormGata ~ Description, data = allData, pch = 16, cex = 0.7, col = c(2:5), add = TRUE);

wilcox.test(subset(allData, allData$Treatment == 1)$NormGata ~ subset(allData, allData$Treatment == 1)$Inner)
wilcox.test(subset(allData, allData$Treatment == 2)$NormGata ~ subset(allData, allData$Treatment == 2)$Inner)

wilcox.test(subset(allData, allData$Inner == FALSE)$NormGata ~ subset(allData, allData$Inner == FALSE)$Treatment)
wilcox.test(subset(allData, allData$Inner == TRUE)$NormGata ~ subset(allData, allData$Inner == TRUE)$Treatment)

dev.off();

pdf(paste("plots", "control_treated_gata3.pdf", sep=.Platform$file.sep), pointsize=13);

boxplot(NormYap ~ Description, data = allData, log = "y", outline = FALSE, ylim = c(0.01,5.0), main = "", ylim = c(0.01,2.0), ylab = "Nuclear YAP1/DAPI Expression", xlab = "", names = c("Outer_Control", "Outer_Treated", "Inner_Control", "Inner_Treated"));
beeswarm(NormYap ~ Description, data = allData, pch = 16, cex = 0.7, col = c(2:4), add = TRUE);

wilcox.test(subset(allData, allData$Inner == FALSE)$NormYap ~ subset(allData, allData$Inner == FALSE)$Treatment)
wilcox.test(subset(allData, allData$Inner == TRUE)$NormYap ~ subset(allData, allData$Inner == TRUE)$Treatment)

dev.off();