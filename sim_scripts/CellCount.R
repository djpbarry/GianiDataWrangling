source("sim_scripts/Definitions.R");

pdf(paste("plots", "sim_centroid_error_versus_snr.pdf", sep=.Platform$file.sep), pointsize=20);
boxplot(Centroid_Error~snr, data=allData, main="", xlab="SNR", ylab="Centroid Error (Microns)", notch=TRUE, col=(c("blue")), ylim=c(0,0.6));

dev.off();

cellCountErrors <- data.frame();

snrs <- unique(allData[[SNR]]);

for(i in snrs){
  temp<-allData[allData$snr == i, ];
  temp<-temp[!duplicated(temp$Index),];
  cellCountErrors <- rbind(cellCountErrors, data.frame(SNR=i,
                                                       Mean_Cell_Count_Error=mean(temp[, CELL_COUNT_ERROR]),
                                                       SD_Cell_Count_Error=sd(temp[, CELL_COUNT_ERROR])));
}

write.csv(cellCountErrors, paste("outputs", "sim_cell_count_errors.csv", sep=.Platform$file.sep));

cellCountErrors <- data.frame();

for(i in snrs){
  temp<-allData[allData$snr == i, ];
  temp<-temp[!duplicated(temp$Index),];
  
  tempData <- data.frame(SNR = factor(rep(i, nrow(temp))),
                         PROP_CELL_COUNT_ERROR = temp$Proportional_Cell_Count_Error);
  
  cellCountErrors <- rbind(cellCountErrors, data.frame(SNR = factor(rep(i, nrow(temp))),
                                                       PROP_CELL_COUNT_ERROR = temp$Proportional_Cell_Count_Error));
}

pdf(paste("plots", "sim_cell_count_errors.pdf", sep=.Platform$file.sep));

p<-ggplot(cellCountErrors, aes(x = PROP_CELL_COUNT_ERROR, fill = SNR)) + geom_histogram(aes(y=0.02*..density..), binwidth = 0.02, alpha =1, position = "dodge");
p<-p + ylab("Relative Frequency");
p<-p + xlab("Cell Count Error");  
p<-p + theme(
  axis.line = element_line(colour = 'black', size = 1.5),
  axis.text.x = element_text(face = "bold", color = "#000000"), axis.text.y = element_text(face = "bold", color = "#000000"),
  legend.position = c(0.1, 0.8),
  text = element_text(size=20,face="bold", colour="black"),
  axis.ticks = element_line(colour = "black", size = 1.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank());
#p<-p + scale_x_discrete(limits=c(0.0,1.0));
p

dev.off();
