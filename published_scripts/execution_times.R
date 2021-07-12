source("published_scripts/definitions.R");
library(ggplot2);

start <- "Start";
input <- "Input:";
proc <- "Processors:";
finish <- "Done";
ncells <- "ncells";

inputDir <- "E:/Dropbox (The Francis Crick)/Debugging/Giani/speed_Test/cells";

files <- list.files(inputDir);

data <- data.frame(matrix(data=0,nrow=length(files),ncol=3));
colnames(data) <- c("Time","Cells", "CPUs");

for(i in 1:length(files)){

slurm <- file(file.path(inputDir,files[i]), "r");

while(TRUE){
  line <- readLines(slurm, n = 1);
  if(length(line) < 1) break;
  lineparts <- strsplit(line, split = " ");
  
  if(grepl(start, lineparts[[1]][1])){
#    print(line);
    startTime <- paste(lineparts[[1]][2],lineparts[[1]][3]);
  } else if (grepl(finish, lineparts[[1]][1])){
 #   print(line);
    endTime <- paste(lineparts[[1]][2],lineparts[[1]][3]);
  } else if (grepl(input, lineparts[[1]][1])){
  #  print(line);
    nCellsStartPos <- regexpr(ncells, lineparts[[1]][3]) + nchar(ncells);
    nCellsEndPos <- nCellsStartPos + regexpr("_", substr(lineparts[[1]][3], nCellsStartPos, nchar(lineparts[[1]][3]))) - 2;
    nCells <- as.numeric(substr(lineparts[[1]][3], nCellsStartPos, nCellsEndPos));
  } else if(grepl(proc, lineparts[[1]][3])){
    ncpus <- as.numeric(lineparts[[1]][4]);
  }
}
if(exists("endTime")){
  timeDiff <- as.numeric(difftime(strptime(endTime, format = "%d/%m/%Y %H:%M:%S"), strptime(startTime, format = "%d/%m/%Y %H:%M:%S"), units="secs"));
  
  data$Time[i] <- timeDiff;
  data$Cells[i] <- nCells;
  data$CPUs[i] <- ncpus;
}
close(slurm);
}

data <- data[data$Time>0,];

ggplot(data, aes(CPUs, Time)) + 
  geom_smooth() +
  geom_point() +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel, plot.title = axislabel) +
  guides(fill = "none", color="none") + xlab("Number of CPUs") + ylab("Execution Time (s)");

data <- data[data$Time<3000,];

ggplot(data, aes(Cells, Time)) + 
  geom_smooth() +
  geom_point() +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.text = axislabel, axis.text.y = axislabel, axis.text.x = axislabel, axis.title.x = axislabel, axis.title.y = axislabel, legend.title = axislabel, plot.title = axislabel) +
  guides(fill = "none", color="none") + xlab("Number of Cells") + ylab("Execution Time (s)");

nrow(data[data$CPUs==16,])
