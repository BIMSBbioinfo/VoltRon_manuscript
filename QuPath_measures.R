# get data
channel_names <- c("CD20", "CD68", "CD21", "FOXP3", "CD45RB")
raw_data <- read.table("data/Tonsil/IF_QuPath/measurements.txt", header = T, sep = "\t")
datax <- raw_data[,grepl("Cell", colnames(raw_data))]
datax <- datax[,grepl("Median", colnames(datax))]
datax <- datax[,grepl(paste(channel_names, collapse = "|"), colnames(datax))]
datax <- datax[,-5]
colnames(datax) <- channel_names

# adjust measures
minmax <- list(
  CD20 = c(284,690),
  CD68 = c(162,374),
  CD21 = c(354,721),
  FOXP3 = c(140,234),
  CD45RB = c(250,479)
)
for(i in 1:length(minmax)){
  datax[,i] <- ifelse(datax[,i] < minmax[[i]][2], datax[,i],minmax[[i]][2])
  datax[,i] <- ifelse(datax[,i] > minmax[[i]][1], datax[,i],minmax[[i]][1])
  datax[,i] <- datax[,i] - min(datax[,i])
}

write.table(datax, file = "data/Tonsil/IF_QuPath/measurements_adjusted.txt")