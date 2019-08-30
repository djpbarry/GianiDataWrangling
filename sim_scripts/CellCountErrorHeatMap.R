source("sim_scripts/Definitions.R");

snrs <- unique(allData[[SNR]]);
cellCounts <- unique(allData[[GROUND_TRUTH_N]]);

snrs<-sort(snrs);
cellCounts<-sort(cellCounts);

cellCountErrors <- data.frame(matrix(data=0.0,nrow=length(cellCounts),ncol=length(snrs)));
colnames(cellCountErrors) <- snrs;
rownames(cellCountErrors) <- cellCounts;

entryCounts <- data.frame(cellCountErrors);

for(i in snrs){
  temp<-allData[allData[[SNR]] == i, ];
  temp<-temp[!duplicated(temp[[INDEX]]),];
  snrCol <- match(i, snrs);
  for(j in 1:nrow(temp)){
    nCells<-temp[[GROUND_TRUTH_N]][j];
    cellRow <- match(nCells, cellCounts);
    cellCountErrors[cellRow, snrCol] <- cellCountErrors[cellRow, snrCol] + temp[[PROP_CELL_COUNT_ERROR]][j];
    entryCounts[cellRow, snrCol] <- entryCounts[cellRow, snrCol] + 1;
  }
}

for(i in 1:nrow(cellCountErrors)){
  for(j in 1:ncol(cellCountErrors)){
    if(entryCounts[i, j] > 0){
      cellCountErrors[i, j] <- cellCountErrors[i, j] / entryCounts[i, j];
    }
  }
}

pdf(paste("plots", "sim_cell_count_errors.pdf", sep=.Platform$file.sep));

#heatmap(as.matrix(cellCountErrors), Colv = NA, Rowv = NA, col=colorRampPalette(brewer.pal(8, "RdYlBu"))(16), xlab="SNR", ylab="Cell Number", cexCol=1.0);

dev.off();

heatMapData <- allData[!duplicated(allData[[INDEX]]),];

ggplot(heatMapData, aes(snr, Ground_Truth_Cell_Count, fill=Proportional_Cell_Count_Error)) + 
  geom_tile() + scale_fill_gradient(low="blue", high="red") + theme_ipsum();
