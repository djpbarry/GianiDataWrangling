controlData <- subset(allData, allData$Treatment == CONTROL_VALUE);
treatedData <- subset(allData, allData$Treatment == TREATED_VALUE);

controlData <- cbind(controlData, data.frame(Rounded_Distance_To_Centre = round(controlData$Distance_To_Centre)));

boxplot(Nucleus_GATA3_Mean_Intensity ~ Rounded_Distance_To_Centre, data = controlData, log = "y", outline = FALSE, ylim = c(2,150));
beeswarm(Nucleus_GATA3_Mean_Intensity ~ Rounded_Distance_To_Centre, data = controlData, pch = 16, cex = 0.7, col = c(2:4), add = TRUE);

boxplot(Nucleus_GATA3_Mean_Intensity ~ Rounded_Distance_To_Centre, data = treatedData, log = "y", outline = FALSE, ylim = c(2,150));
beeswarm(Nucleus_GATA3_Mean_Intensity ~ Rounded_Distance_To_Centre, data = treatedData, pch = 16, cex = 0.7, col = c(2:4), add = TRUE);