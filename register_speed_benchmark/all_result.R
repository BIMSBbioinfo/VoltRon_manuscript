library(RBioFormats)
library(peakRAM)
library(bench)
library(magick)
library(dplyr)
library(ggplot2)
library(scales)

####
# functions ####
####

convert_to_bytes <- function(x) {
  
  # Remove spaces for safety
  x <- trimws(x)
  
  if (grepl("MB$", x, ignore.case = TRUE)) {
    value <- as.numeric(sub("MB$", "", x, ignore.case = TRUE))
    return(value * 1024^2)  # MB → bytes
  } else if (grepl("GB$", x, ignore.case = TRUE)) {
    value <- as.numeric(sub("GB$", "", x, ignore.case = TRUE))
    return(value * 1024^3)  # GB → bytes
  } else if (grepl("KB$", x, ignore.case = TRUE)) {
    value <- as.numeric(sub("KB$", "", x, ignore.case = TRUE))
    return(value * 1024)  # KB → bytes
  } 
  
  return(x)
}

transform_breaks_time <- function(x){
  breaks <- c(-1,c(0.5,1,30,70,120))
  x_ind <- as.integer(cut(x, breaks = breaks))
  diff <- breaks[-1] - breaks[-length(breaks)]
  low <- breaks[x_ind]
  high <- breaks[x_ind + 1]
  x_norm <- ((x - low)/(high - low)) + (x_ind - 1)
  return(x_norm)
}

transform_breaks_mem <- function(x){
  breaks <- c(1,1024, 1024^2, 1024^3, 4*1024^3)
  x_ind <- as.integer(cut(x, breaks = breaks))
  diff <- breaks[-1] - breaks[-length(breaks)]
  low <- breaks[x_ind]
  high <- breaks[x_ind + 1]
  # x_norm <- ((x - low)/(high - low)) + (x_ind - 1)
  x_norm <- (x)/(1024^(x_ind)) + (x_ind - 1)
  return(x_norm)
}

####
# import ####
####

# get results
list_files <- list.files("results/", full.names = TRUE)
list_files_names <- basename(list_files)
list_files_names <- gsub("\\.[^.]*$", "", list_files_names)
ind <- c(1,2,4,3,5)
list_files <- list_files[ind]
list_files_names <- list_files_names[ind]
all_results <- NULL
for(i in 1:length(list_files_names)){
  cur_res <- read.csv(list_files[i])
  colnames(cur_res) <- c("measure", list_files_names[i])
  if(i > 1){
    all_results <- all_results %>% full_join(cur_res)
  } else {
    all_results <- cur_res
  }
}
all_results

####
# adjust visualize data ####
####

# adjust data for visualization
all_results$data <- rep(c("image", "observations", "image", "molecules"), 
                        c(6,3,1,3))
# all_results <- reshape2::melt(all_results, id.vars = c("data", "measure"))
# all_results <- na.omit(all_results)

####
# data types ####
####

## image ####
all_results_image <- all_results[all_results$data == "image",]
rownames(all_results_image) <- all_results_image$measure
all_results_image <- t(all_results_image[,c(-1, -1*ncol(all_results_image))])
all_results_image <- as.data.frame(all_results_image)
all_results_image$name <- paste0(rownames(all_results_image), " \n (", 
                                 all_results_image$size_ref, ") vs (", all_results_image$size_query, ")")
all_results_image <- all_results_image[,c(-1,-2)]
all_results_image <- reshape2::melt(all_results_image, id.vars = "name")
all_results_image$value <- sapply(all_results_image$value, convert_to_bytes)

## cells ####
all_results_cells <- all_results[all_results$data == "observations",]
rownames(all_results_cells) <- all_results_cells$measure
all_results_cells <- t(all_results_cells[,c(-1, -1*ncol(all_results_cells))])
all_results_cells <- as.data.frame(all_results_cells)
all_results_cells$name <- paste0(rownames(all_results_cells), " (", all_results_cells$n_cells, ")")
all_results_cells <- all_results_cells[,-1]
all_results_cells <- reshape2::melt(all_results_cells, id.vars = "name")
all_results_cells$value <- sapply(all_results_cells$value, convert_to_bytes)

## molecules ####
all_results_mol <- all_results[all_results$data == "molecules",]
rownames(all_results_mol) <- all_results_mol$measure
all_results_mol <- t(all_results_mol[,c(-1, -1*ncol(all_results_mol))])
all_results_mol <- as.data.frame(all_results_mol)
all_results_mol$name <- paste0(rownames(all_results_mol), " (", all_results_mol$n_molecules, ")")
all_results_mol <- all_results_mol[,-1]
all_results_mol <- reshape2::melt(all_results_mol, id.vars = "name")
all_results_mol$value <- sapply(all_results_mol$value, convert_to_bytes)
all_results_mol <- na.omit(all_results_mol)

## concatenate ####
all_results <- rbind(
  data.frame(all_results_image, type = "image"),
  data.frame(all_results_cells, type = "cells/spots"),
  data.frame(all_results_mol, type = "molecules")
)
tmp <- all_results$type
tmp <- factor(tmp, levels = c("cells/spots", "molecules", "image"))
all_results$type <- tmp
all_results$measure <- sapply(as.character(all_results$variable), 
                           \(x) strsplit(x, split = '_')[[1]][1])
all_results$value <- as.numeric(all_results$value)
tmp <- all_results$name
lbls <- levels(factor(all_results$name))
lbls <- lbls[c(2,4,9,6,12,10,7,1,3,8,5,11)]
tmp <- factor(tmp, levels = lbls)
all_results$name <- tmp

####
# visualize time ####
####

# visualize time 
all_results_time <- all_results[all_results$measure %in% "time",]
all_results_time$name <- droplevels(all_results_time$name)
breaks <- c(0.5,1,30,70,120)
all_results_time$value <- transform_breaks_time(all_results_time$value)
all_results_time$data <- sapply(as.character(all_results_time$name), \(x){
  strsplit(x, split = " ")[[1]][1]
})
ggplot(mapping = aes(x = name, y = value, fill = data), data = all_results_time) + 
  geom_bar(stat = "identity") + 
  facet_grid(~type, scales = "free_x", space = "free_x") +
  scale_y_continuous(
    limits = c(0, length(breaks)),
    breaks = 0:length(breaks),   # <-- choose your custom break points
    labels = as.character(c(0,breaks))
  ) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(legend.position = "top") + 
  ylab("Time (seconds)") + xlab("")
  
####
# visualize memory ####
####

all_results_mem <- all_results[!all_results$measure %in% "time",]
all_results_mem$name <- droplevels(all_results_mem$name)
breaks <- c(1,1024, 1024^2, 1024^3, 4*1024^3)
# all_results_mem$value <- log(all_results_mem$value)/log(1024)
# all_results_mem$value2 <- transform_breaks_mem(all_results_mem$value)
all_results_mem$data <- sapply(as.character(all_results_mem$name), \(x){
  strsplit(x, split = " ")[[1]][1]
})
all_results_mem$type <- ifelse(all_results_mem$type != "image", 
                               as.character(all_results_mem$type), 
                               ifelse(all_results_mem$variable == "mem_reg",
                                      "Image (Reg.)",
                                      ifelse(all_results_mem$variable == "mem_warp",
                                             "Image (Warp)", 
                                             "Image (Size)")))
tmp <- all_results_mem$type
tmp <- factor(tmp, levels = c("cells/spots", "molecules", "Image (Reg.)", "Image (Warp)", "Image (Size)"))
all_results_mem$type <- tmp
ggplot(mapping = aes(x = name, y = value, fill = data), data = all_results_mem) + 
  geom_bar(stat = "identity") + 
  facet_grid(~type, scales = "free_x", space = "free_x") +
  scale_y_continuous(
    trans = "log2",
    breaks = c(1, 1024, 1024^2, 1024^3, 4*1024^3),
    labels = c("1B", "1KB", "1MB", "1GB", "4GB"),
    limits = c(1, 4*1024^3)
  ) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(legend.position = "top") + 
  ylab("Memory Consumption (Bytes)") + xlab("")