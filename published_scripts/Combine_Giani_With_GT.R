source("published_scripts/definitions.R");
library(ggplot2);

groundTruthData <- read.csv('D:/Dropbox (The Francis Crick)/Debugging/Giani/Ground_Truth.csv');
gianiData <- read.csv('D:/Dropbox (The Francis Crick)/Debugging/Giani/GIANI-3.1.2_Output_with_Ground_Truth.csv');
imarisData <- read.csv('D:/Dropbox (The Francis Crick)/Debugging/Giani/Sim_Image_snr1.000000_ncells2000_Statistics/Sim_Image_snr1.000000_ncells2000_Cell_Position_with_Ground_Truth.csv');
gianiColNames <- colnames(gianiData);
nucGianiData <-  subset(gianiData, grepl(NUCLEUS, gianiData[[LABEL]]));
nGD <- nrow(nucGianiData);
nGT <- nrow(groundTruthData);
appendage <- data.frame(matrix(
  data = NaN,
  nrow = nGT,
  ncol = 8
))
colnames(appendage) <-
  c(
    M_NUC_CENTROID_X,
    M_NUC_CENTROID_Y,
    M_NUC_CENTROID_Z,
    IM_NUC_CENTROID_X,
    IM_NUC_CENTROID_Y,
    IM_NUC_CENTROID_Z,
    GIANI_NUC_CENTROID_ERROR,
    IMARIS_NUC_CENTROID_ERROR
  )

groundTruthData$Nearest_Neighbour <- .Machine$integer.max;

for(j in 1:nGT){
  gtx1 <- groundTruthData[[GT_CENTROID_X]][j];
  gty1 <- groundTruthData[[GT_CENTROID_Y]][j];
  gtz1 <- groundTruthData[[GT_CENTROID_Z]][j];
  for(i in 1:nGT){
    if(!(i == j)){
      gtx2 <- groundTruthData[[GT_CENTROID_X]][i];
      gty2 <- groundTruthData[[GT_CENTROID_Y]][i];
      gtz2 <- groundTruthData[[GT_CENTROID_Z]][i];
      dist <- euclidDistance(c(gtx1, gty1, gtz1), c(gtx2, gty2, gtz2));
      if(dist < groundTruthData$Nearest_Neighbour[j]){
        groundTruthData$Nearest_Neighbour[j] <- dist;
        groundTruthData$Nearest_Neighbour_Index[j] <- groundTruthData$Ground_Truth_Cell_Index[i];
      }
    }
  }
}

allData <- cbind(groundTruthData, appendage);
errorData <- data.frame();

detections <- c('Neither', 'Giani Only', 'Imaris Only', 'Both');

allData$Detected <- detections[1];

for (i in 1:nGT) {
  gtIndex <- allData$Ground_Truth_Cell_Index[i];
  giani_Index <- match(gtIndex, nucGianiData$Ground_Truth_Cell_Index);
  imaris_Index <- match(gtIndex, imarisData$Ground_Truth_Cell_Index);
  
  gtx <- allData[[GT_CENTROID_X]][i];
  gty <- allData[[GT_CENTROID_Y]][i];
  gtz <- allData[[GT_CENTROID_Z]][i];
  
  if(is.finite(giani_Index)){
    giani_x <- nucGianiData[[GD_CENTROID_X]][giani_Index];
    giani_y <- nucGianiData[[GD_CENTROID_Y]][giani_Index];
    giani_z <- nucGianiData[[GD_CENTROID_Z]][giani_Index];
      
    dist <- euclidDistance(c(giani_x, giani_y, giani_z), c(gtx, gty, gtz));
      
    allData[[M_NUC_CENTROID_X]][gtIndex] <- giani_x;
    allData[[M_NUC_CENTROID_Y]][gtIndex] <- giani_y;
    allData[[M_NUC_CENTROID_Z]][gtIndex] <- giani_z;
    allData[[GIANI_NUC_CENTROID_ERROR]][gtIndex] <- dist;
    errorData <- rbind(errorData, data.frame(type='GIANI', error=dist, Nearest_Neighbour=allData$Nearest_Neighbour[i]));
    allData$Detected[i] <- detections[2];
  }
  
  if(is.finite(imaris_Index)){
    imaris_x <- imarisData[[IMARIS_X]][imaris_Index];
    imaris_y <- imarisData[[IMARIS_Y]][imaris_Index];
    imaris_z <- imarisData[[IMARIS_Z]][imaris_Index];
    
    dist <- euclidDistance(c(imaris_x, imaris_y, imaris_z), c(gtx, gty, gtz));
    
    allData[[IM_NUC_CENTROID_X]][gtIndex] <- imaris_x;
    allData[[IM_NUC_CENTROID_Y]][gtIndex] <- imaris_y;
    allData[[IM_NUC_CENTROID_Z]][gtIndex] <- imaris_z;
    allData[[IMARIS_NUC_CENTROID_ERROR]][gtIndex] <- dist;
    errorData <- rbind(errorData, data.frame(type='IMARIS', error=dist, Nearest_Neighbour=allData$Nearest_Neighbour[i]));
    if(allData$Detected[i] == detections[2]){
      allData$Detected[i] <- detections[4];
    } else{
      allData$Detected[i] <- detections[3];
    }
  }
  
}

filteredGianiData <- allData[!is.nan(allData$GIANI_Nuc_Centroid_Error),];
filteredImarisData <- allData[!is.nan(allData$Imaris_Nuc_Centroid_Error),];

allData$Giani_Detected <- is.finite(allData$GIANI_Nuc_Centroid_Error);
allData$Imaris_Detected <- is.finite(allData$Imaris_Nuc_Centroid_Error);

ns <- c(sum(allData$Detected==detections[4]),
        sum(allData$Detected==detections[2]),
        sum(allData$Detected==detections[3]),
        sum(allData$Detected==detections[1]));

quantile(allData[allData$Detected=='Neither',]$Nearest_Neighbour, c(.25,.5,.75,.95))
quantile(allData[allData$Detected=='Imaris Only',]$Nearest_Neighbour, c(.25,.5,.75))
quantile(allData[allData$Detected=='Giani Only',]$Nearest_Neighbour, c(.25,.5,.75))

ggplot(allData, aes(x=Detected, y=Nearest_Neighbour, color=Detected)) + 
#  geom_violin(scale="count",bw=1,size=1,width=3,alpha=0.5) +
  geom_text(x=1,y=95, label=paste("n = ", ns[1]),color=black,size=10) +
  geom_text(x=2,y=55,label=paste("n = ", ns[2]),color=black,size=10) +
  geom_text(x=3,y=20,label=paste("n = ", ns[3]),color=black,size=10) +
  geom_text(x=4,y=20,label=paste("n = ", ns[4]),color=black,size=10) +
  geom_jitter(width=0.15,alpha=0.5,size=3) +
  geom_boxplot(width=0.6,color=black,coef=20,alpha=0.0) +
  scale_color_manual(values=c(cc[2],cc[1],cc[4],cc[3])) +
  ylim(c(0,100)) +
#p <- p + geom_dotplot(binaxis="y", stackdir="center", binwidth = 0.6);
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel, plot.title = axislabel) +
  guides(fill = "none", color="none") + xlab(element_blank()) + ylab("Distance to Nearest\nNeighbour (\U03BCm)") +
  scale_x_discrete(breaks=unique(allData$Detected), labels = c("Both", "GIANI\nOnly", "Imaris\nOnly", "Neither"));

ggplot(errorData, aes(x=error, fill=type)) +
  geom_density(alpha=.35) +
  #geom_histogram( color=black, alpha=0.5, position = 'identity', binwidth = 0.2) +
  scale_fill_manual(values=c(blue,orange)) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel, plot.title = axislabel) +
  theme(legend.position = c(0.8,0.8)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  labs(fill="") +
  xlab("Localisation Error (\U03BCm)") + ylab("Normalised Frequency") +
  xlim(c(0,7));

#ggplot(errorData, aes(x=type, y=error, fill=type)) +
#  geom_boxplot(width=0.1,color=black,coef=20) +
#  geom_violin(bw=0.1,size=1,alpha=0.5) +
#  scale_fill_manual(values=c(blue,orange)) +
#  ylim(c(0,7)) +
#  theme_linedraw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#  theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel, plot.title = axislabel) +
#  guides(fill = "none") + xlab("") + ylab("Nucleus Localisation Error (\U03BCm)");

ggplot(errorData, aes(x=error, y=Nearest_Neighbour, color=type)) + 
  geom_point(size=2,alpha=0.35) +
  scale_color_manual(values=c(blue,orange)) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(),legend.position=c(0.85,0.2),legend.background = element_blank(), legend.box.background = element_rect(colour = "black"), legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, plot.title = axislabel) +
#  theme(legend.position = c(0.9,0.15)) +
  scale_y_log10("Distance to Nearest\nNeighbour (\U03BCm)") +
  scale_x_log10(limits=c(0.05,200)) +
  labs(color="") +
  xlab("Localisation Error (\U03BCm)");

filename <- "all_advanced_giani_data.csv";

if(simpleGiani){
  filename <- "all_simple_giani_data.csv";
}

write.csv(allData, file.path("outputs", filename));
