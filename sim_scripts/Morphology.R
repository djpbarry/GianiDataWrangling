source("sim_scripts/Definitions.R");

pdf(paste("plots", "volume_error_versus_snr.pdf", sep=.Platform$file.sep));
boxplot(Volume_Error~snr, data=allData, main="Volume Error Versus SNR", xlab="SNR", ylab="Volume Error (microns^3)", notch=TRUE, col=(c("darkgreen")));
legend("topright", c("2"), fill=c("darkgreen"));

dev.off();

pdf(paste("plots", "norm_volume_error_versus_snr.pdf", sep=.Platform$file.sep));
summary<-boxplot(Normalised_Volume_Error~snr, data=allData, main="Normalised Volume Error Versus SNR", xlab="SNR", ylab="Normalised Volume Error", notch=TRUE, plot=FALSE, col=(c("darkgreen")));
legend("topright", c("2"), fill=c("darkgreen"));

dev.off();

pdf(paste("plots", "volume_error_versus_cell_number.pdf", sep=.Platform$file.sep));
boxplot(Volume_Error~Ground_Truth_Cell_Count, data=allData, main="Volume Error Versus Cell Count", xlab="Cell Count", ylab="Volume Error (microns^3)", notch=TRUE, col=(c("darkgreen")));
legend("topright", c("2"), fill=c("darkgreen"));

dev.off();

pdf(paste("plots", "norm_volume_error_versus_cell_number.pdf", sep=.Platform$file.sep));
boxplot(Normalised_Volume_Error~Ground_Truth_Cell_Count, data=allData, main="Normalised Volume Error Versus Cell Count", xlab="Cell Count", ylab="Normalised Volume Error", notch=TRUE, col=(c("darkgreen")));
legend("topright", c("2"), fill=c("darkgreen"));

dev.off();

pdf(paste("plots", "norm_volume_error_versus_z", sep=.Platform$file.sep));

plot(allData[allData$snr==1,names(allData) %in% c(GT_CENTROID_Z)], allData[allData$snr==1,names(allData) %in% c(NORM_VOL_ERROR)], xlab="Z (microns)", main="Normalised Volume Error Versus Z", ylab="Normalised Volume Error", ylim=c(0,1), xlim=c(20,120), col="red", pch=15);
points(allData[allData$snr==2,names(allData) %in% c(GT_CENTROID_Z)], allData[allData$snr==2,names(allData) %in% c(NORM_VOL_ERROR)], col="green", pch=15);
points(allData[allData$snr==3,names(allData) %in% c(GT_CENTROID_Z)], allData[allData$snr==3,names(allData) %in% c(NORM_VOL_ERROR)], col="blue", pch=15);
legend("topright", c("1","2","3"), fill=c("red", "green", "blue"));

dev.off();
