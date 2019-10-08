dist_thresh = 0.75;

controlData <- subset(allData, allData$Treatment == CONTROL_VALUE);
treatedData <- subset(allData, allData$Treatment == TREATED_VALUE);

controlData$Inner <- controlData$Distance_To_Centre < dist_thresh;
treatedData$Inner <- treatedData$Distance_To_Centre < dist_thresh;

controlData <- cbind(controlData, data.frame(Norm_gata = controlData$Nucleus_GATA3_Mean_Intensity / controlData$Nucleus_DAPI_Mean_Intensity));
treatedData <- cbind(treatedData, data.frame(Norm_gata = treatedData$Nucleus_GATA3_Mean_Intensity / treatedData$Nucleus_DAPI_Mean_Intensity));

controlData <- cbind(controlData, data.frame(Norm_yap = controlData$Nucleus_YAP1_Mean_Intensity / controlData$Nucleus_DAPI_Mean_Intensity));
treatedData <- cbind(treatedData, data.frame(Norm_yap = treatedData$Nucleus_YAP1_Mean_Intensity / treatedData$Nucleus_DAPI_Mean_Intensity));

boxplot(Norm_gata ~ Inner, data = controlData, log = "y", outline = FALSE, main = "norm_gata_control", ylim = c(0.01,2.0));
beeswarm(Norm_gata ~ Inner, data = controlData, pch = 16, cex = 0.7, col = c(2:4), add = TRUE);

boxplot(Norm_gata ~ Inner, data = treatedData, log = "y", outline = FALSE, main = "norm_gata_treated", ylim = c(0.01,2.0));
beeswarm(Norm_gata ~ Inner, data = treatedData, pch = 16, cex = 0.7, col = c(2:4), add = TRUE);

boxplot(Norm_yap ~ Inner, data = controlData, log = "y", outline = FALSE, main = "norm_yap_control", ylim = c(0.01,2.0));
beeswarm(Norm_yap ~ Inner, data = controlData, pch = 16, cex = 0.7, col = c(2:4), add = TRUE);

boxplot(Norm_yap ~ Inner, data = treatedData, log = "y", outline = FALSE, main = "norm_yap_treated", ylim = c(0.01,2.0));
beeswarm(Norm_yap ~ Inner, data = treatedData, pch = 16, cex = 0.7, col = c(2:4), add = TRUE);
