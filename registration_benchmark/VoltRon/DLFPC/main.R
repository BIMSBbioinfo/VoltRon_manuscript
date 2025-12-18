# library(VoltRon)
library(patchwork)
library(anndataR)
library(anndata)
source("../../allresults/DLFPC/accuracy.R")

####
# import ####
####

# import datasets
DLPFC_1 <- importVisium("../../../../../data/10X_Visium_DLPFC/151673/", sample_name = "DLPFC_1")
DLPFC_2 <- importVisium("../../../../../data/10X_Visium_DLPFC/151674", sample_name = "DLPFC_2")
DLPFC_3 <- importVisium("../../../../../data/10X_Visium_DLPFC/151675", sample_name = "DLPFC_3")
DLPFC_4 <- importVisium("../../../../../data/10X_Visium_DLPFC/151676", sample_name = "DLPFC_4")
DLPFC_507 <- importVisium("../../../../../data/10X_Visium_DLPFC/151507/", sample_name = "DLPFC_507")
DLPFC_508 <- importVisium("../../../../../data/10X_Visium_DLPFC/151508/", sample_name = "DLPFC_508")
DLPFC_509 <- importVisium("../../../../../data/10X_Visium_DLPFC/151509/", sample_name = "DLPFC_509")
DLPFC_510 <- importVisium("../../../../../data/10X_Visium_DLPFC/151510/", sample_name = "DLPFC_510")
DLPFC_669 <- importVisium("../../../../../data/10X_Visium_DLPFC/151669/", sample_name = "DLPFC_669")
DLPFC_670 <- importVisium("../../../../../data/10X_Visium_DLPFC/151670/", sample_name = "DLPFC_670")
DLPFC_671 <- importVisium("../../../../../data/10X_Visium_DLPFC/151671/", sample_name = "DLPFC_671")
DLPFC_672 <- importVisium("../../../../../data/10X_Visium_DLPFC/151672/", sample_name = "DLPFC_672")

# get annotation
DLPFC_1$layer_guess <- as.character(readRDS("../../allresults/DLFPC/151673_labels.rds")[gsub("_Assay1","", vrSpatialPoints(DLPFC_1))])
DLPFC_2$layer_guess <- as.character(readRDS("../../allresults/DLFPC/151674_labels.rds")[gsub("_Assay1","", vrSpatialPoints(DLPFC_2))])
DLPFC_3$layer_guess <- as.character(readRDS("../../allresults/DLFPC/151675_labels.rds")[gsub("_Assay1","", vrSpatialPoints(DLPFC_3))])
DLPFC_4$layer_guess <- as.character(readRDS("../../allresults/DLFPC/151676_labels.rds")[gsub("_Assay1","", vrSpatialPoints(DLPFC_4))])
DLPFC_507$layer_guess <- as.character(readRDS("../../allresults/DLFPC/151507_labels.rds")[gsub("_Assay1","", vrSpatialPoints(DLPFC_507))])
DLPFC_508$layer_guess <- as.character(readRDS("../../allresults/DLFPC/151508_labels.rds")[gsub("_Assay1","", vrSpatialPoints(DLPFC_508))])
DLPFC_509$layer_guess <- as.character(readRDS("../../allresults/DLFPC/151509_labels.rds")[gsub("_Assay1","", vrSpatialPoints(DLPFC_509))])
DLPFC_510$layer_guess <- as.character(readRDS("../../allresults/DLFPC/151510_labels.rds")[gsub("_Assay1","", vrSpatialPoints(DLPFC_510))])
DLPFC_669$layer_guess <- as.character(readRDS("../../allresults/DLFPC/151669_labels.rds")[gsub("_Assay1","", vrSpatialPoints(DLPFC_669))])
DLPFC_670$layer_guess <- as.character(readRDS("../../allresults/DLFPC/151670_labels.rds")[gsub("_Assay1","", vrSpatialPoints(DLPFC_670))])
DLPFC_671$layer_guess <- as.character(readRDS("../../allresults/DLFPC/151671_labels.rds")[gsub("_Assay1","", vrSpatialPoints(DLPFC_671))])
DLPFC_672$layer_guess <- as.character(readRDS("../../allresults/DLFPC/151672_labels.rds")[gsub("_Assay1","", vrSpatialPoints(DLPFC_672))])

####
# align ####
####

# # 1000 20
# DLPFC_list <- list(DLPFC_1, DLPFC_2)
# # reg1and2 <- registerSpatialData(object_list = DLPFC_list)
# # saveRDS(reg1and2$mapping_parameters, file = "parameters/reg1and2_parameters.rds")
# reg1and2 <- registerSpatialData(object_list = DLPFC_list, mapping_parameters = readRDS(file = "parameters/reg1and2_parameters.rds"), interactive = FALSE)
# # 800 50
# DLPFC_list <- list(DLPFC_3, DLPFC_4)
# # reg3and4 <- registerSpatialData(object_list = DLPFC_list)
# # saveRDS(reg3and4$mapping_parameters, file = "parameters/reg3and4_parameters.rds")
# reg3and4 <- registerSpatialData(object_list = DLPFC_list, mapping_parameters = readRDS(file = "parameters/reg3and4_parameters.rds"), interactive = FALSE)
# merge_list <- c(reg1and2$registered_spat, reg3and4$registered_spat)
# SRBlock <- merge(merge_list[[1]], merge_list[-1], samples = "DLPFC_Block")
# SRBlock

####
# align with crops ####
####

####
## crop data ####
####

# DLPFC_1_cropped <- subset(DLPFC_1, interactive = TRUE)
# saveRDS(DLPFC_1_cropped$subset_info_list[[1]], file = "parameters/DLPFC_1_cropinfo.rds")
DLPFC_1_cropped <- subset(DLPFC_1, image = readRDS("parameters/DLPFC_1_cropinfo.rds"))
rm(DLPFC_1)

# DLPFC_2_cropped <- subset(DLPFC_2, interactive = TRUE)
# saveRDS(DLPFC_2_cropped$subset_info_list[[1]], file = "parameters/DLPFC_2_cropinfo.rds")
DLPFC_2_cropped <- subset(DLPFC_2, image = readRDS("parameters/DLPFC_2_cropinfo.rds"))
rm(DLPFC_2)

# DLPFC_3_cropped <- subset(DLPFC_3, interactive = TRUE)
# saveRDS(DLPFC_3_cropped$subset_info_list[[1]], file = "parameters/DLPFC_3_cropinfo.rds")
DLPFC_3_cropped <- subset(DLPFC_3, image = readRDS("parameters/DLPFC_3_cropinfo.rds"))
rm(DLPFC_3)

# DLPFC_4_cropped <- subset(DLPFC_4, interactive = TRUE)
# saveRDS(DLPFC_4_cropped$subset_info_list[[1]], file = "parameters/DLPFC_4_cropinfo.rds")
DLPFC_4_cropped <- subset(DLPFC_4, image = readRDS("parameters/DLPFC_4_cropinfo.rds"))
rm(DLPFC_4)

# DLPFC_507_cropped <- subset(DLPFC_507, interactive = TRUE)
# saveRDS(DLPFC_507_cropped$subset_info_list[[1]], file = "parameters/DLPFC_507_cropinfo.rds")
DLPFC_507_cropped <- subset(DLPFC_507, image = readRDS("parameters/DLPFC_507_cropinfo.rds"))
rm(DLPFC_507)

# DLPFC_508_cropped <- subset(DLPFC_508, interactive = TRUE)
# saveRDS(DLPFC_508_cropped$subset_info_list[[1]], file = "parameters/DLPFC_508_cropinfo.rds")
DLPFC_508_cropped <- subset(DLPFC_508, image = readRDS("parameters/DLPFC_508_cropinfo.rds"))
rm(DLPFC_508)

# DLPFC_509_cropped <- subset(DLPFC_509, interactive = TRUE)
# saveRDS(DLPFC_509_cropped$subset_info_list[[1]], file = "parameters/DLPFC_509_cropinfo.rds")
DLPFC_509_cropped <- subset(DLPFC_509, image = readRDS("parameters/DLPFC_509_cropinfo.rds"))
# rm(DLPFC_509)

# DLPFC_510_cropped <- subset(DLPFC_510, interactive = TRUE)
# saveRDS(DLPFC_510_cropped$subset_info_list[[1]], file = "parameters/DLPFC_510_cropinfo.rds")
DLPFC_510_cropped <- subset(DLPFC_510, image = readRDS("parameters/DLPFC_510_cropinfo.rds"))
# rm(DLPFC_510)

# DLPFC_669_cropped <- subset(DLPFC_669, interactive = TRUE)
# saveRDS(DLPFC_669_cropped$subset_info_list[[1]], file = "parameters/DLPFC_669_cropinfo.rds")
DLPFC_669_cropped <- subset(DLPFC_669, image = readRDS("parameters/DLPFC_669_cropinfo.rds"))
rm(DLPFC_669)

# DLPFC_670_cropped <- subset(DLPFC_670, interactive = TRUE)
# saveRDS(DLPFC_670_cropped$subset_info_list[[1]], file = "parameters/DLPFC_670_cropinfo.rds")
DLPFC_670_cropped <- subset(DLPFC_670, image = readRDS("parameters/DLPFC_670_cropinfo.rds"))
rm(DLPFC_670)

# DLPFC_671_cropped <- subset(DLPFC_671, interactive = TRUE)
# saveRDS(DLPFC_671_cropped$subset_info_list[[1]], file = "parameters/DLPFC_671_cropinfo.rds")
DLPFC_671_cropped <- subset(DLPFC_671, image = readRDS("parameters/DLPFC_671_cropinfo.rds"))
rm(DLPFC_671)

# DLPFC_672_cropped <- subset(DLPFC_672, interactive = TRUE)
# saveRDS(DLPFC_672_cropped$subset_info_list[[1]], file = "parameters/DLPFC_672_cropinfo.rds")
DLPFC_672_cropped <- subset(DLPFC_672, image = readRDS("parameters/DLPFC_672_cropinfo.rds"))
rm(DLPFC_672)

####
## align ####
####

# 1500 20
DLPFC_list <- list(DLPFC_1_cropped, DLPFC_2_cropped)
# reg1and2 <- registerSpatialData(object_list = DLPFC_list)
# saveRDS(reg1and2$mapping_parameters, file = "parameters/reg1and2_parameters_cropped.rds")
mapping_parameters = readRDS(file = "parameters/reg1and2_parameters_cropped.rds")
# mapping_parameters$MAX_FEATURES <- as.character(1500)
reg1and2 <- registerSpatialData(object_list = DLPFC_list, mapping_parameters = mapping_parameters, interactive = FALSE)
# reg1and2 <- registerSpatialData(object_list = DLPFC_list, mapping_parameters = mapping_parameters)

# 800 50
DLPFC_list <- list(DLPFC_3_cropped, DLPFC_4_cropped)
# reg3and4 <- registerSpatialData(object_list = DLPFC_list)
# saveRDS(reg3and4$mapping_parameters, file = "parameters/reg3and4_parameters_cropped.rds")
mapping_parameters = readRDS(file = "parameters/reg3and4_parameters_cropped.rds")
reg3and4 <- registerSpatialData(object_list = DLPFC_list, mapping_parameters = mapping_parameters, interactive = FALSE)

# merge
merge_list <- c(reg1and2$registered_spat, reg3and4$registered_spat)
SRBlock <- merge(merge_list[[1]], merge_list[-1], samples = "DLPFC_Block")
SRBlock

# homography
DLPFC_list <- list(DLPFC_507_cropped, DLPFC_508_cropped)
# reg7and8 <- registerSpatialData(object_list = DLPFC_list)
# saveRDS(reg7and8$mapping_parameters, file = "parameters/reg7and8_parameters_cropped.rds")
mapping_parameters = readRDS(file = "parameters/reg7and8_parameters_cropped.rds")
reg7and8 <- registerSpatialData(object_list = DLPFC_list, mapping_parameters = mapping_parameters, interactive = FALSE)

# homography
# DLPFC_list <- list(DLPFC_509_cropped, DLPFC_510_cropped)
DLPFC_list <- list(DLPFC_509, DLPFC_510)
# reg9and10 <- registerSpatialData(object_list = DLPFC_list)
# saveRDS(reg9and10$mapping_parameters, file = "parameters/reg9and10_parameters_cropped.rds")
mapping_parameters = readRDS(file = "parameters/reg9and10_parameters_cropped.rds")
reg9and10 <- registerSpatialData(object_list = DLPFC_list, mapping_parameters = mapping_parameters, interactive = FALSE)

# merge
merge_list <- c(reg7and8$registered_spat, reg9and10$registered_spat)
SRBlock2 <- merge(merge_list[[1]], merge_list[-1], samples = "DLPFC_Block")
SRBlock2

# homography
DLPFC_list <- list(DLPFC_669_cropped, DLPFC_670_cropped)
# reg69and70 <- registerSpatialData(object_list = DLPFC_list)
# saveRDS(reg69and70$mapping_parameters, file = "parameters/reg69and70_parameters_cropped.rds")
mapping_parameters = readRDS(file = "parameters/reg69and70_parameters_cropped.rds")
reg69and70 <- registerSpatialData(object_list = DLPFC_list, mapping_parameters = mapping_parameters, interactive = FALSE)
# reg69and70 <- registerSpatialData(object_list = DLPFC_list, mapping_parameters = mapping_parameters)

# homography
DLPFC_list <- list(DLPFC_671_cropped, DLPFC_672_cropped)
# reg71and72 <- registerSpatialData(object_list = DLPFC_list)
# saveRDS(reg71and72$mapping_parameters, file = "parameters/reg71and72_parameters_cropped.rds")
mapping_parameters = readRDS(file = "parameters/reg71and72_parameters_cropped.rds")
reg71and72 <- registerSpatialData(object_list = DLPFC_list, mapping_parameters = mapping_parameters, interactive = FALSE)

# merge
merge_list <- c(reg69and70$registered_spat, reg71and72$registered_spat)
SRBlock3 <- merge(merge_list[[1]], merge_list[-1], samples = "DLPFC_Block")
SRBlock3

####
# alignment accuracy ####
####

# 1 and 2
obj1 <- vrCoordinates(SRBlock, assay = "Assay1")[,c("x", "y")]
obj2 <- vrCoordinates(SRBlock, assay = "Assay2")[,c("x", "y")]
obj2_orig <- vrCoordinates(SRBlock, assay = "Assay2", spatial_name = "main")[,c("x", "y")]
spots <- vrSpatialPoints(SRBlock, assay = "Assay1")
label1 <- Metadata(SRBlock, assay = "Assay1")[['layer_guess']]
label1[is.na(label1)] <- "unlabeled"
label2 <- Metadata(SRBlock, assay = "Assay2")[['layer_guess']]
label2[is.na(label2)] <- "unlabeled"
dist <- vrAssayParams(SRBlock[["Assay1"]])[['nearestpost.distance']]
results <- get_adjacent_label(obj1, obj2, label1, label2, dist)
result_orig <- get_adjacent_label(obj1, obj2_orig, label1, label2, dist)
saveRDS(data.frame(label = label1, results = results, result_orig = result_orig, row.names = spots), file = "results/1and2/predictions.rds")

# 3 and 4
obj1 <- vrCoordinates(SRBlock, assay = "Assay3")[,c("x", "y")]
obj2 <- vrCoordinates(SRBlock, assay = "Assay4")[,c("x", "y")]
obj2_orig <- vrCoordinates(SRBlock, assay = "Assay4", spatial_name = "main")[,c("x", "y")]
spots <- vrSpatialPoints(SRBlock, assay = "Assay3")
label1 <- Metadata(SRBlock, assay = "Assay3")[['layer_guess']]
label1[is.na(label1)] <- "unlabeled"
label2 <- Metadata(SRBlock, assay = "Assay4")[['layer_guess']]
label2[is.na(label2)] <- "unlabeled"
dist <- vrAssayParams(SRBlock[["Assay3"]])[['nearestpost.distance']]
results <- get_adjacent_label(obj1, obj2, label1, label2, dist)
result_orig <- get_adjacent_label(obj1, obj2_orig, label1, label2, dist)
saveRDS(data.frame(label = label1, results = results, result_orig = result_orig, row.names = spots), file = "results/3and4/predictions.rds")

# 7 and 8
obj1 <- vrCoordinates(SRBlock2, assay = "Assay1", reg = TRUE)[,c("x", "y")]
obj2 <- vrCoordinates(SRBlock2, assay = "Assay2", reg = TRUE)[,c("x", "y")]
obj2_orig <- vrCoordinates(SRBlock2, assay = "Assay2", spatial_name = "main")[,c("x", "y")]
spots <- vrSpatialPoints(SRBlock2, assay = "Assay1")
label1 <- Metadata(SRBlock2, assay = "Assay1")[['layer_guess']]
label1[is.na(label1)] <- "unlabeled"
label2 <- Metadata(SRBlock2, assay = "Assay2")[['layer_guess']]
label2[is.na(label2)] <- "unlabeled"
dist <- vrAssayParams(SRBlock2[["Assay1"]])[['nearestpost.distance']]
results <- get_adjacent_label(obj1, obj2, label1, label2, dist)
result_orig <- get_adjacent_label(obj1, obj2_orig, label1, label2, dist)
saveRDS(data.frame(label = label1, results = results, result_orig = result_orig, row.names = spots), file = "results/7and8/predictions.rds")

# 9 and 10
obj1 <- vrCoordinates(SRBlock2, assay = "Assay3", reg = TRUE)[,c("x", "y")]
obj2 <- vrCoordinates(SRBlock2, assay = "Assay4", reg = TRUE)[,c("x", "y")]
obj2_orig <- vrCoordinates(SRBlock2, assay = "Assay4", spatial_name = "main")[,c("x", "y")]
spots <- vrSpatialPoints(SRBlock2, assay = "Assay3")
label1 <- Metadata(SRBlock2, assay = "Assay3")[['layer_guess']]
label1[is.na(label1)] <- "unlabeled"
label2 <- Metadata(SRBlock2, assay = "Assay4")[['layer_guess']]
label2[is.na(label2)] <- "unlabeled"
dist <- vrAssayParams(SRBlock2[["Assay3"]])[['nearestpost.distance']]
results <- get_adjacent_label(obj1, obj2, label1, label2, dist)
result_orig <- get_adjacent_label(obj1, obj2_orig, label1, label2, dist)
saveRDS(data.frame(label = label1, results = results, result_orig = result_orig, row.names = spots), file = "results/9and10/predictions.rds")

# 69 and 70
obj1 <- vrCoordinates(SRBlock3, assay = "Assay1", reg = TRUE)[,c("x", "y")]
obj2 <- vrCoordinates(SRBlock3, assay = "Assay2", reg = TRUE)[,c("x", "y")]
obj2_orig <- vrCoordinates(SRBlock3, assay = "Assay2", spatial_name = "main")[,c("x", "y")]
spots <- vrSpatialPoints(SRBlock3, assay = "Assay1")
label1 <- Metadata(SRBlock3, assay = "Assay1")[['layer_guess']]
label1[is.na(label1)] <- "unlabeled"
label2 <- Metadata(SRBlock3, assay = "Assay2")[['layer_guess']]
label2[is.na(label2)] <- "unlabeled"
dist <- vrAssayParams(SRBlock3[["Assay1"]])[['nearestpost.distance']]
results <- get_adjacent_label(obj1, obj2, label1, label2, dist)
result_orig <- get_adjacent_label(obj1, obj2_orig, label1, label2, dist)
saveRDS(data.frame(label = label1, results = results, result_orig = result_orig, row.names = spots), file = "results/69and70/predictions.rds")

# 71 and 72
obj1 <- vrCoordinates(SRBlock3, assay = "Assay3", reg = TRUE)[,c("x", "y")]
obj2 <- vrCoordinates(SRBlock3, assay = "Assay4", reg = TRUE)[,c("x", "y")]
obj2_orig <- vrCoordinates(SRBlock3, assay = "Assay4", spatial_name = "main")[,c("x", "y")]
spots <- vrSpatialPoints(SRBlock3, assay = "Assay3")
label1 <- Metadata(SRBlock3, assay = "Assay3")[['layer_guess']]
label1[is.na(label1)] <- "unlabeled"
label2 <- Metadata(SRBlock3, assay = "Assay4")[['layer_guess']]
label2[is.na(label2)] <- "unlabeled"
dist <- vrAssayParams(SRBlock3[["Assay3"]])[['nearestpost.distance']]
results <- get_adjacent_label(obj1, obj2, label1, label2, dist)
result_orig <- get_adjacent_label(obj1, obj2_orig, label1, label2, dist)
saveRDS(data.frame(label = label1, results = results, result_orig = result_orig, row.names = spots), file = "results/71and72/predictions.rds")

####
# visualize ####
####

####
## 1 vs 2 ####
####

set.seed(1)
td1 <- tempfile(fileext = ".png")
colors <- setNames(object = hue_pal(8), nm = unique(SRBlock$layer_guess))
g1 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay1", colors = colors)
ggsave(filename = td1, plot = g1, width = 8, height = 8)
td2 <- tempfile(fileext = ".png")
g2 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay2", spatial = "main_reg", colors = colors)
ggsave(filename = td2, plot = g2, width = 8, height = 8)
c(image_read(td1), image_read(td2)) %>% 
  magick::image_join()
c(rep(image_read(td1),5), rep(image_read(td2),5)) %>% 
  magick::image_join() %>% magick::image_write_gif(path = "temp.gif")

set.seed(1)
td1 <- tempfile(fileext = ".png")
colors <- setNames(object = hue_pal(8), nm = unique(SRBlock$layer_guess))
g1 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay1", colors = colors)
ggsave(filename = td1, plot = g1, width = 8, height = 8)
td2 <- tempfile(fileext = ".png")
g2 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay2", spatial = "main", colors = colors)
ggsave(filename = td2, plot = g2, width = 8, height = 8)
c(image_read(td1), image_read(td2)) %>% 
  magick::image_join()

####
## 3 vs 4 ####
####

set.seed(1)
td1 <- tempfile(fileext = ".png")
colors <- setNames(object = hue_pal(8), nm = unique(SRBlock$layer_guess))
g1 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay3", colors = colors)
ggsave(filename = td1, plot = g1, width = 8, height = 8)
td2 <- tempfile(fileext = ".png")
g2 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay4", spatial = "main_reg", colors = colors)
ggsave(filename = td2, plot = g2, width = 8, height = 8)
c(image_read(td1), image_read(td2)) %>% 
  magick::image_join()

set.seed(1)
td1 <- tempfile(fileext = ".png")
colors <- setNames(object = hue_pal(8), nm = unique(SRBlock$layer_guess))
g1 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay3", colors = colors)
ggsave(filename = td1, plot = g1, width = 8, height = 8)
td2 <- tempfile(fileext = ".png")
g2 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay4", spatial = "main", colors = colors)
ggsave(filename = td2, plot = g2, width = 8, height = 8)
c(image_read(td1), image_read(td2)) %>% 
  magick::image_join()

####
## 7 vs 8 ####
####

set.seed(1)
td1 <- tempfile(fileext = ".png")
colors <- setNames(object = hue_pal(8), nm = unique(SRBlock2$layer_guess))
g1 <- vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay1", colors = colors)
ggsave(filename = td1, plot = g1, width = 8, height = 8)
td2 <- tempfile(fileext = ".png")
g2 <- vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay2", spatial = "main_reg", colors = colors)
ggsave(filename = td2, plot = g2, width = 8, height = 8)
c(image_read(td1), image_read(td2)) %>% 
  magick::image_join()

####
## 69 vs 70 ####
####

set.seed(1)
td1 <- tempfile(fileext = ".png")
colors <- setNames(object = hue_pal(6), nm = unique(SRBlock3$layer_guess))
g1 <- vrSpatialPlot(SRBlock3, group.by = "layer_guess", assay = "Assay1", colors = colors)
ggsave(filename = td1, plot = g1, width = 8, height = 8)
td2 <- tempfile(fileext = ".png")
g2 <- vrSpatialPlot(SRBlock3, group.by = "layer_guess", assay = "Assay2", spatial = "main_reg", colors = colors)
ggsave(filename = td2, plot = g2, width = 8, height = 8)
c(image_read(td1), image_read(td2)) %>% 
  magick::image_join()

####
## 71 vs 72 ####
####

set.seed(1)
td1 <- tempfile(fileext = ".png")
colors <- setNames(object = hue_pal(6), nm = unique(SRBlock3$layer_guess))
g1 <- vrSpatialPlot(SRBlock3, group.by = "layer_guess", assay = "Assay3", colors = colors)
ggsave(filename = td1, plot = g1, width = 8, height = 8)
td2 <- tempfile(fileext = ".png")
g2 <- vrSpatialPlot(SRBlock3, group.by = "layer_guess", assay = "Assay4", spatial = "main_reg", colors = colors)
ggsave(filename = td2, plot = g2, width = 8, height = 8)
c(image_read(td1), image_read(td2)) %>% 
  magick::image_join()

set.seed(1)
td1 <- tempfile(fileext = ".png")
colors <- setNames(object = hue_pal(6), nm = unique(SRBlock3$layer_guess))
g1 <- vrSpatialPlot(SRBlock3, group.by = "layer_guess", assay = "Assay3", colors = colors)
ggsave(filename = td1, plot = g1, width = 8, height = 8)
td2 <- tempfile(fileext = ".png")
g2 <- vrSpatialPlot(SRBlock3, group.by = "layer_guess", assay = "Assay4", spatial = "main", colors = colors)
ggsave(filename = td2, plot = g2, width = 8, height = 8)
c(image_read(td1), image_read(td2)) %>% 
  magick::image_join()

####
## visualize all ####
####

layer_guess <- c("Layer3", "Layer1", "WM", "Layer5", "Layer6", "Layer2", "Layer4", NA)
colors <- setNames(object = hue_pal(8), nm = layer_guess)
colors <- as.list(colors)
g1 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay1", colors = colors) |>
  addSpatialLayer(SRBlock, group.by = "layer_guess", assay = "Assay2", alpha = 0.5, spatial = "main", colors = colors) + 
  theme(legend.position = "none") + labs(title = "")
g2 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay1", colors = colors) |>
  addSpatialLayer(SRBlock, group.by = "layer_guess", assay = "Assay2", alpha = 0.5, spatial = "main_reg", colors = colors) +
  theme(legend.position = "none") + labs(title = "")
g3 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay3", colors = colors) |>
  addSpatialLayer(SRBlock, group.by = "layer_guess", assay = "Assay4", alpha = 0.5, spatial = "main", colors = colors) +
  theme(legend.position = "none") + labs(title = "")
g4 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay3", colors = colors) |>
  addSpatialLayer(SRBlock, group.by = "layer_guess", assay = "Assay4", alpha = 0.5, spatial = "main_reg", colors = colors) +
  theme(legend.position = "none") + labs(title = "")
g5 <- vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay1", colors = colors) |>
  addSpatialLayer(SRBlock2, group.by = "layer_guess", assay = "Assay2", alpha = 0.5, spatial = "main", colors = colors) + 
  theme(legend.position = "none") + labs(title = "")
g6 <- vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay1", colors = colors) |>
  addSpatialLayer(SRBlock2, group.by = "layer_guess", assay = "Assay2", alpha = 0.5, spatial = "main_reg", colors = colors) + 
  theme(legend.position = "none") + labs(title = "")
g7 <- vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay3", crop = TRUE, colors = colors) |>
  addSpatialLayer(SRBlock2, group.by = "layer_guess", assay = "Assay4", alpha = 0.5, spatial = "main", colors = colors) +
  theme(legend.position = "none") + labs(title = "")
g8 <- vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay3", crop = TRUE, colors = colors) |>
  addSpatialLayer(SRBlock2, group.by = "layer_guess", assay = "Assay4", alpha = 0.5, spatial = "main_reg", colors = colors) +
  theme(legend.position = "none") + labs(title = "")
colors <- colors[unique(SRBlock3$layer_guess)]
g9 <- vrSpatialPlot(SRBlock3, group.by = "layer_guess", assay = "Assay1", colors = colors) |>
  addSpatialLayer(SRBlock3, group.by = "layer_guess", assay = "Assay2", alpha = 0.5, spatial = "main", colors = colors) +
  theme(legend.position = "none") + labs(title = "")
g10 <- vrSpatialPlot(SRBlock3, group.by = "layer_guess", assay = "Assay1", colors = colors) |>
  addSpatialLayer(SRBlock3, group.by = "layer_guess", assay = "Assay2", alpha = 0.5, spatial = "main_reg", colors = colors) + 
  theme(legend.position = "none") + labs(title = "")
g11 <- vrSpatialPlot(SRBlock3, group.by = "layer_guess", assay = "Assay3", colors = colors) |>
  addSpatialLayer(SRBlock3, group.by = "layer_guess", assay = "Assay4", alpha = 0.5, spatial = "main", colors = colors) +
  theme(legend.position = "none") + labs(title = "")
g12 <- vrSpatialPlot(SRBlock3, group.by = "layer_guess", assay = "Assay3", colors = colors) |>
  addSpatialLayer(SRBlock3, group.by = "layer_guess", assay = "Assay4", alpha = 0.5, spatial = "main_reg", colors = colors) + 
  theme(legend.position = "none") + labs(title = "")

ggpubr::ggarrange(plotlist = list(g1,g3,g5,g7,g9,g11,g2,g4,g6,g8,g10,g12), nrow = 2, ncol = 6)
ggsave("results/voltron_alloverlay.pdf", width = 36, height = 13)

vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay1", colors = colors) + theme(legend.position = "top")
ggsave("results/voltron_labelled.pdf", width = 6, height = 7)
