source("sim_scripts/Definitions.R");

allData <- read.csv("E:/Dropbox (The Francis Crick)/Giani_Imaris_Comp/CompiledGianiData.csv")

red <- "#DD0000";
green <- "#008800";
yellow <- "#FFFF00";

snrs <- unique(allData[[SNR]]);
cellCounts <- unique(allData[[GROUND_TRUTH_N]]);

snrs<-sort(snrs);
cellCounts<-sort(cellCounts);

cellCountErrors <- data.frame();

for(i in snrs){
  temp1<-allData[allData[[SNR]] == i, ];
  for(j in cellCounts){
    temp2 <- temp1[temp1$Ground_Truth_Cell_Count == j, ];
    runs <- unique(temp2$run_);
    errors <- matrix(0,3,3);
    for(r in runs){
      temp3 <- temp2[temp2$run_ == r, ];
      #temp4 <- temp3[is.finite(temp3$Centroid_Error),];
      temp4 <- temp3[is.finite(temp3$Nuc_Centroid_Error),];
      errors[r,1] <- mean(temp4$Cell_Count_Error);
      errors[r,2] <- mean(temp4$Proportional_Cell_Count_Error);
      #errors[r,3] <- mean(temp4$Centroid_Error);
      errors[r,3] <- mean(temp4$Nuc_Centroid_Error);
    }
    #errors <- temp2[is.finite(temp2$Centroid_Error),];
    #cellCountErrors <- rbind(cellCountErrors, c(i, j, mean(abs(errors$Cell_Count_Error)), mean(abs(errors$Proportional_Cell_Count_Error)), mean(errors$Nuc_Centroid_Error), mean(errors$Cell_Centroid_Error), mean(abs(errors$Normalised_Volume_Error)), mean(abs(errors$Proportionate_Volume_Error))));
    cellCountErrors <- rbind(cellCountErrors, c(i, j, mean(abs(errors[,1])), mean(abs(errors[,2])), mean(errors[,3])));
  }
}

colnames(cellCountErrors) <- c(SNR, "N_CELLS", "CELL_COUNT_ERROR", "PROP_CELL_COUNT_ERROR", "MEAN_NUC_CENTROID_ERROR");

#pdf(paste("plots", "sim_cell_count_errors.pdf", sep=.Platform$file.sep));

#heatmap(as.matrix(cellCountErrors), Colv = NA, Rowv = NA, col=colorRampPalette(brewer.pal(8, "RdYlBu"))(16), xlab="SNR", ylab="Cell Number", cexCol=1.0);

#dev.off();

#heatMapData <- allData[!duplicated(allData[[INDEX]]),];
axislabel <- element_text(hjust=0.5, size=25, colour = "black");

#Cell count plot
p <- ggplot(cellCountErrors, aes(as.character(snr), N_CELLS));
p <- p + geom_raster(aes(fill = CELL_COUNT_ERROR),interpolate=FALSE);
p <- p + scale_fill_gradient2(low=green, high=red, mid = yellow, midpoint = 2.5, name = expression(paste("E" [c])), limits = c(0,3), breaks =c(0,1,2));
p <- p + xlab("SNR") + ylab("Number of Cells") + labs(fill = "Cell Count Error");
p <- p + scale_x_discrete(labels=c("0.10", "0.25", "0.50", "0.75", "1.00", "2.00", "3.00"));
p <- p + theme_minimal();
p + theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel);

#Centroid error plot
p <- ggplot(cellCountErrors, aes(as.character(snr), N_CELLS));
p <- p + geom_raster(aes(fill = MEAN_NUC_CENTROID_ERROR),interpolate=FALSE);
p <- p + scale_fill_gradient2(low=green, high=red, mid = yellow, midpoint = 0.5, name = expression(paste("E" [nl])), limits = c(0, 1), breaks =c(0, 0.5, 1));
p <- p + xlab("SNR") + ylab("Number of Cells") + labs(fill = "Localisation Error");
p <- p + scale_x_discrete(labels=c("0.10", "0.25", "0.50", "0.75", "1.00", "2.00", "3.00"));
p <- p + theme_minimal();
p + theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel);

#Cell centroid error
p <- ggplot(cellCountErrors, aes(as.character(snr), N_CELLS));
p <- p + geom_raster(aes(fill = MEAN_CELL_CENTROID_ERROR),interpolate=FALSE);
p <- p + scale_fill_gradient2(low=green, high=red, mid = yellow, midpoint = 6, name = expression(paste("E" [cl])), limits = c(1, 11), breaks =c(0, 5, 10));
p <- p + xlab("SNR") + ylab("Number of Cells") + labs(fill = "Localisation Error");
p <- p + scale_x_discrete(labels=c("0.10", "0.25", "0.50", "0.75", "1.00", "2.00", "3.00"));
p <- p + theme_minimal();
p + theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel);





#Proportional cell count plot
p <- ggplot(cellCountErrors, aes(as.character(snr), N_CELLS));
p <- p + geom_raster(aes(fill = PROP_CELL_COUNT_ERROR),interpolate=FALSE);
p <- p + scale_fill_gradient2(low=green, high=red, mid = yellow, midpoint = 0.03, name = expression(paste("E" [cc])), limits = c(0.00, 0.06), breaks =c(0.00, 0.02, 0.04, 0.06));
p <- p + xlab("SNR") + ylab("Number of Cells") + labs(fill = "Cell Count Error");
p <- p + scale_x_discrete(labels=c("0.10", "0.25", "0.50", "0.75", "1.00", "2.00", "3.00"));
p <- p + theme_minimal();
p + theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel);

#Segmentation error
p <- ggplot(cellCountErrors, aes(as.character(snr), N_CELLS));
p <- p + geom_raster(aes(fill = VOLUME_ERROR),interpolate=FALSE);
p <- p + scale_fill_gradient2(low=green, high=red, mid = yellow, midpoint = 0.295, name = expression(paste("E" [cs])), limits = c(0.12, 0.47), breaks =c(0.1, 0.2, 0.3, 0.4));
p <- p + xlab("SNR") + ylab("Number of Cells") + labs(fill = "Ecs");
p <- p + scale_x_discrete(labels=c("0.10", "0.25", "0.50", "0.75", "1.00", "2.00", "3.00"));
p <- p + theme_minimal();
p + theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel);

#Corrected segmentation error
p <- ggplot(cellCountErrors, aes(as.character(snr), N_CELLS));
p <- p + geom_raster(aes(fill = CORRECTED_VOLUME_ERROR),interpolate=FALSE);
p <- p + scale_fill_gradient2(low=green, high=red, mid = yellow, midpoint=0.012, name = expression(paste("E" [ms])), limits = c(0.002, 0.022), breaks =c(0.005, 0.010, 0.015, 0.020));
p <- p + xlab("SNR") + ylab("Number of Cells") + labs(fill = "Ems");
p <- p + scale_x_discrete(labels=c("0.10", "0.25", "0.50", "0.75", "1.00", "2.00", "3.00"));
p <- p + theme_minimal();
p + theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel);

#plot(x, controlData, xlab=xLabel, main=title, ylab=yLabel, ylim=yLimits, col="red", pch=15);
plot(allData$Normalised_Distance_to_Centre, allData$Normalised_Volume_Error);

mean(allData$Centroid_Error, na.rm = TRUE)
