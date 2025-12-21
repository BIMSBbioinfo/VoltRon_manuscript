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
# get aligned ####
####

####
## 1,2,3,4 ####
####

# temp <- anndata::read_h5ad("results/151673adata_aligned.h5ad")
# temp2 <- anndata::read_h5ad("data/151673.h5ad")

temp <- anndata::read_h5ad("results/151673adata_aligned.h5ad")
coords <- temp$obsm$spatial
rownames(coords) <- sapply(temp$obs_names, \(x) strsplit(x, split = "\\.")[[1]][1])
rownames(coords) <- paste0(rownames(coords), "_", vrAssayNames(DLPFC_1, assay = "Assay1"))
DLPFC_1 <- subset(DLPFC_1, spatialpoints = rownames(coords))
vrImages(DLPFC_1[["Assay1"]], reg = TRUE) <- formImage(coords = coords)

temp <- anndata::read_h5ad("results/151674adata_aligned.h5ad")
coords <- temp$obsm$spatial
rownames(coords) <- sapply(temp$obs_names, \(x) strsplit(x, split = "\\.")[[1]][1])
rownames(coords) <- paste0(rownames(coords), "_", vrAssayNames(DLPFC_2, assay = "Assay1"))
DLPFC_2 <- subset(DLPFC_2, spatialpoints = rownames(coords))
vrImages(DLPFC_2[["Assay1"]], reg = TRUE) <- formImage(coords = coords)

temp <- anndata::read_h5ad("results/151675adata_aligned.h5ad")
coords <- temp$obsm$spatial
rownames(coords) <- sapply(temp$obs_names, \(x) strsplit(x, split = "\\.")[[1]][1])
rownames(coords) <- paste0(rownames(coords), "_", vrAssayNames(DLPFC_3, assay = "Assay1"))
DLPFC_3 <- subset(DLPFC_3, spatialpoints = rownames(coords))
vrImages(DLPFC_3[["Assay1"]], reg = TRUE) <- formImage(coords = coords)

temp <- anndata::read_h5ad("results/151676adata_aligned.h5ad")
coords <- temp$obsm$spatial
rownames(coords) <- sapply(temp$obs_names, \(x) strsplit(x, split = "\\.")[[1]][1])
rownames(coords) <- paste0(rownames(coords), "_", vrAssayNames(DLPFC_4, assay = "Assay1"))
DLPFC_4 <- subset(DLPFC_4, spatialpoints = rownames(coords))
vrImages(DLPFC_4[["Assay1"]], reg = TRUE) <- formImage(coords = coords)

merge_list <- list(DLPFC_1, DLPFC_2, DLPFC_3, DLPFC_4)
rm(DLPFC_1, DLPFC_2, DLPFC_3, DLPFC_4)
SRBlock <- merge(merge_list[[1]], merge_list[-1], samples = "DLPFC_Block")
SRBlock

####
## 7,8,9,10 ####
####

temp <- anndata::read_h5ad("results/151507adata_aligned.h5ad")
coords <- temp$obsm$spatial
rownames(coords) <- sapply(temp$obs_names, \(x) strsplit(x, split = "\\.")[[1]][1])
rownames(coords) <- paste0(rownames(coords), "_", vrAssayNames(DLPFC_507, assay = "Assay1"))
DLPFC_507 <- subset(DLPFC_507, spatialpoints = rownames(coords))
vrImages(DLPFC_507[["Assay1"]], reg = TRUE) <- formImage(coords = coords)

temp <- anndata::read_h5ad("results/151508adata_aligned.h5ad")
coords <- temp$obsm$spatial
rownames(coords) <- sapply(temp$obs_names, \(x) strsplit(x, split = "\\.")[[1]][1])
rownames(coords) <- paste0(rownames(coords), "_", vrAssayNames(DLPFC_508, assay = "Assay1"))
DLPFC_508 <- subset(DLPFC_508, spatialpoints = rownames(coords))
vrImages(DLPFC_508[["Assay1"]], reg = TRUE) <- formImage(coords = coords)

temp <- anndata::read_h5ad("results/151509adata_aligned.h5ad")
coords <- temp$obsm$spatial
rownames(coords) <- sapply(temp$obs_names, \(x) strsplit(x, split = "\\.")[[1]][1])
rownames(coords) <- paste0(rownames(coords), "_", vrAssayNames(DLPFC_509, assay = "Assay1"))
DLPFC_509 <- subset(DLPFC_509, spatialpoints = rownames(coords))
vrImages(DLPFC_509[["Assay1"]], reg = TRUE) <- formImage(coords = coords)

temp <- anndata::read_h5ad("results/151510adata_aligned.h5ad")
coords <- temp$obsm$spatial
rownames(coords) <- sapply(temp$obs_names, \(x) strsplit(x, split = "\\.")[[1]][1])
rownames(coords) <- paste0(rownames(coords), "_", vrAssayNames(DLPFC_510, assay = "Assay1"))
DLPFC_510 <- subset(DLPFC_510, spatialpoints = rownames(coords))
vrImages(DLPFC_510[["Assay1"]], reg = TRUE) <- formImage(coords = coords)

merge_list <- list(DLPFC_507, DLPFC_508, DLPFC_509, DLPFC_510)
rm(DLPFC_507, DLPFC_508, DLPFC_509, DLPFC_510)
SRBlock2 <- merge(merge_list[[1]], merge_list[-1], samples = "DLPFC_Block")
SRBlock2

####
## 69,70,71,72 ####
####

temp <- anndata::read_h5ad("results/151669adata_aligned.h5ad")
coords <- temp$obsm$spatial
rownames(coords) <- sapply(temp$obs_names, \(x) strsplit(x, split = "\\.")[[1]][1])
rownames(coords) <- paste0(rownames(coords), "_", vrAssayNames(DLPFC_669, assay = "Assay1"))
DLPFC_669 <- subset(DLPFC_669, spatialpoints = rownames(coords))
vrImages(DLPFC_669[["Assay1"]], reg = TRUE) <- formImage(coords = coords)

temp <- anndata::read_h5ad("results/151670adata_aligned.h5ad")
coords <- temp$obsm$spatial
rownames(coords) <- sapply(temp$obs_names, \(x) strsplit(x, split = "\\.")[[1]][1])
rownames(coords) <- paste0(rownames(coords), "_", vrAssayNames(DLPFC_670, assay = "Assay1"))
DLPFC_670 <- subset(DLPFC_670, spatialpoints = rownames(coords))
vrImages(DLPFC_670[["Assay1"]], reg = TRUE) <- formImage(coords = coords)

temp <- anndata::read_h5ad("results/151671adata_aligned.h5ad")
coords <- temp$obsm$spatial
rownames(coords) <- sapply(temp$obs_names, \(x) strsplit(x, split = "\\.")[[1]][1])
rownames(coords) <- paste0(rownames(coords), "_", vrAssayNames(DLPFC_671, assay = "Assay1"))
DLPFC_671 <- subset(DLPFC_671, spatialpoints = rownames(coords))
vrImages(DLPFC_671[["Assay1"]], reg = TRUE) <- formImage(coords = coords)

temp <- anndata::read_h5ad("results/151672adata_aligned.h5ad")
coords <- temp$obsm$spatial
rownames(coords) <- sapply(temp$obs_names, \(x) strsplit(x, split = "\\.")[[1]][1])
rownames(coords) <- paste0(rownames(coords), "_", vrAssayNames(DLPFC_672, assay = "Assay1"))
DLPFC_672 <- subset(DLPFC_672, spatialpoints = rownames(coords))
vrImages(DLPFC_672[["Assay1"]], reg = TRUE) <- formImage(coords = coords)

merge_list <- list(DLPFC_669, DLPFC_670, DLPFC_671, DLPFC_672)
rm(DLPFC_669, DLPFC_670, DLPFC_671, DLPFC_672)
SRBlock3 <- merge(merge_list[[1]], merge_list[-1], samples = "DLPFC_Block")
SRBlock3

####
# alignment accuracy ####
####

# 1 and 2
obj1 <- vrCoordinates(SRBlock, assay = "Assay1", reg = TRUE)[,c("x", "y")]
obj2 <- vrCoordinates(SRBlock, assay = "Assay2", reg = TRUE)[,c("x", "y")]
spots <- vrSpatialPoints(SRBlock, assay = "Assay1")
label1 <- Metadata(SRBlock, assay = "Assay1")[['layer_guess']]
label1[is.na(label1)] <- "unlabeled"
label2 <- Metadata(SRBlock, assay = "Assay2")[['layer_guess']]
label2[is.na(label2)] <- "unlabeled"
dist <- vrAssayParams(SRBlock[["Assay1"]])[['nearestpost.distance']]
results <- get_adjacent_label(obj1, obj2, label1, label2, dist)
saveRDS(data.frame(label = label1, results = results, row.names = spots), file = "results/1and2/predictions.rds")

# 3 and 4
obj1 <- vrCoordinates(SRBlock, assay = "Assay3", reg = TRUE)[,c("x", "y")]
obj2 <- vrCoordinates(SRBlock, assay = "Assay4", reg = TRUE)[,c("x", "y")]
spots <- vrSpatialPoints(SRBlock, assay = "Assay3")
label1 <- Metadata(SRBlock, assay = "Assay3")[['layer_guess']]
label1[is.na(label1)] <- "unlabeled"
label2 <- Metadata(SRBlock, assay = "Assay4")[['layer_guess']]
label2[is.na(label2)] <- "unlabeled"
dist <- vrAssayParams(SRBlock[["Assay3"]])[['nearestpost.distance']]
results <- get_adjacent_label(obj1, obj2, label1, label2, dist)
saveRDS(data.frame(label = label1, results = results, row.names = spots), file = "results/3and4/predictions.rds")

# 7 and 8
obj1 <- vrCoordinates(SRBlock2, assay = "Assay1", reg = TRUE)[,c("x", "y")]
obj2 <- vrCoordinates(SRBlock2, assay = "Assay2", reg = TRUE)[,c("x", "y")]
spots <- vrSpatialPoints(SRBlock2, assay = "Assay1")
label1 <- Metadata(SRBlock2, assay = "Assay1")[['layer_guess']]
label1[is.na(label1)] <- "unlabeled"
label2 <- Metadata(SRBlock2, assay = "Assay2")[['layer_guess']]
label2[is.na(label2)] <- "unlabeled"
dist <- vrAssayParams(SRBlock2[["Assay1"]])[['nearestpost.distance']]
results <- get_adjacent_label(obj1, obj2, label1, label2, dist)
saveRDS(data.frame(label = label1, results = results, row.names = spots), file = "results/7and8/predictions.rds")

# 9 and 10
obj1 <- vrCoordinates(SRBlock2, assay = "Assay3", reg = TRUE)[,c("x", "y")]
obj2 <- vrCoordinates(SRBlock2, assay = "Assay4", reg = TRUE)[,c("x", "y")]
spots <- vrSpatialPoints(SRBlock2, assay = "Assay3")
label1 <- Metadata(SRBlock2, assay = "Assay3")[['layer_guess']]
label1[is.na(label1)] <- "unlabeled"
label2 <- Metadata(SRBlock2, assay = "Assay4")[['layer_guess']]
label2[is.na(label2)] <- "unlabeled"
dist <- vrAssayParams(SRBlock2[["Assay3"]])[['nearestpost.distance']]
results <- get_adjacent_label(obj1, obj2, label1, label2, dist)
saveRDS(data.frame(label = label1, results = results, row.names = spots), file = "results/9and10/predictions.rds")

# 69 and 70
obj1 <- vrCoordinates(SRBlock3, assay = "Assay1", reg = TRUE)[,c("x", "y")]
obj2 <- vrCoordinates(SRBlock3, assay = "Assay2", reg = TRUE)[,c("x", "y")]
spots <- vrSpatialPoints(SRBlock3, assay = "Assay1")
label1 <- Metadata(SRBlock3, assay = "Assay1")[['layer_guess']]
label1[is.na(label1)] <- "unlabeled"
label2 <- Metadata(SRBlock3, assay = "Assay2")[['layer_guess']]
label2[is.na(label2)] <- "unlabeled"
dist <- vrAssayParams(SRBlock3[["Assay1"]])[['nearestpost.distance']]
results <- get_adjacent_label(obj1, obj2, label1, label2, dist)
saveRDS(data.frame(label = label1, results = results, row.names = spots), file = "results/69and70/predictions.rds")

# 71 and 72
obj1 <- vrCoordinates(SRBlock3, assay = "Assay3", reg = TRUE)[,c("x", "y")]
obj2 <- vrCoordinates(SRBlock3, assay = "Assay4", reg = TRUE)[,c("x", "y")]
spots <- vrSpatialPoints(SRBlock3, assay = "Assay3")
label1 <- Metadata(SRBlock3, assay = "Assay3")[['layer_guess']]
label1[is.na(label1)] <- "unlabeled"
label2 <- Metadata(SRBlock3, assay = "Assay4")[['layer_guess']]
label2[is.na(label2)] <- "unlabeled"
dist <- vrAssayParams(SRBlock3[["Assay3"]])[['nearestpost.distance']]
results <- get_adjacent_label(obj1, obj2, label1, label2, dist)
saveRDS(data.frame(label = label1, results = results, row.names = spots), file = "results/71and72/predictions.rds")

####
# visualize ####
####

####
## 1 vs 2 ####
####

set.seed(1)
td1 <- tempfile(fileext = ".png")
colors <- setNames(object = hue_pal(7), nm = unique(SRBlock$layer_guess))
g1 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay1", spatial = "main_reg", colors = colors)
ggsave(filename = td1, plot = g1, width = 8, height = 8)
td2 <- tempfile(fileext = ".png")
g2 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay1", spatial = "main_reg", colors = colors, alpha = 0) |>
  addSpatialLayer(SRBlock, group.by = "layer_guess", assay = "Assay2", spatial = "main_reg", colors = colors)
ggsave(filename = td2, plot = g2, width = 8, height = 8)
c(image_read(td1), image_read(td2)) %>% 
  magick::image_join() %>% magick::image_write_gif(path = "temp.gif")

# g1 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay1", spatial = "main_reg") |>
#   addSpatialLayer(SRBlock, group.by = "layer_guess", assay = "Assay2", spatial = "main_reg")
vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay1", spatial = "main_reg") |>
  addSpatialLayer(SRBlock, group.by = "layer_guess", assay = "Assay2", spatial = "main_reg", alpha = 0.5)

####
## 3 vs 4 ####
####

set.seed(1)
td1 <- tempfile(fileext = ".png")
colors <- setNames(object = hue_pal(7), nm = unique(SRBlock$layer_guess))
g1 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay3", spatial = "main_reg", colors = colors)
ggsave(filename = td1, plot = g1, width = 8, height = 8)
td2 <- tempfile(fileext = ".png")
g2 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay3", spatial = "main_reg", colors = colors, alpha = 0) |>
  addSpatialLayer(SRBlock, group.by = "layer_guess", assay = "Assay4", spatial = "main_reg", colors = colors)
ggsave(filename = td2, plot = g2, width = 8, height = 8)
c(image_read(td1), image_read(td2)) %>% 
  magick::image_join()

# vrSpatialPlot(SRBlock, group.by = "Sample", assay = "Assay3", spatial = "main_reg") |>
#   addSpatialLayer(SRBlock, group.by = "Sample", assay = "Assay4", spatial = "main_reg", colors = list(DLPFC_Block = "blue"))
vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay3", spatial = "main_reg") |>
  addSpatialLayer(SRBlock, group.by = "layer_guess", assay = "Assay4", spatial = "main_reg", alpha = 0.5)

vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay3", spatial = "main_reg") |
vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay4", spatial = "main_reg")

####
## 7 vs 8 ####
####

set.seed(1)
td1 <- tempfile(fileext = ".png")
colors <- setNames(object = hue_pal(7), nm = unique(SRBlock2$layer_guess))
g1 <- vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay1", spatial = "main_reg", colors = colors)
ggsave(filename = td1, plot = g1, width = 8, height = 8)
td2 <- tempfile(fileext = ".png")
g2 <- vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay1", spatial = "main_reg", colors = colors, alpha = 0) |>
  addSpatialLayer(SRBlock2, group.by = "layer_guess", assay = "Assay2", spatial = "main_reg", colors = colors)
ggsave(filename = td2, plot = g2, width = 8, height = 8)
c(image_read(td1), image_read(td2)) %>% 
  magick::image_join()

vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay1", spatial = "main_reg") |>
  addSpatialLayer(SRBlock2, group.by = "layer_guess", assay = "Assay2", spatial = "main_reg", alpha = 0.5)

####
## 9 vs 10 ####
####

set.seed(1)
td1 <- tempfile(fileext = ".png")
colors <- setNames(object = hue_pal(7), nm = unique(SRBlock2$layer_guess))
g1 <- vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay3", spatial = "main_reg", colors = colors)
ggsave(filename = td1, plot = g1, width = 8, height = 8)
td2 <- tempfile(fileext = ".png")
g2 <- vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay3", spatial = "main_reg", colors = colors, alpha = 0) |>
  addSpatialLayer(SRBlock2, group.by = "layer_guess", assay = "Assay4", spatial = "main_reg", colors = colors)
ggsave(filename = td2, plot = g2, width = 8, height = 8)
c(image_read(td1), image_read(td2)) %>% 
  magick::image_join()

vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay3", spatial = "main_reg") |>
  addSpatialLayer(SRBlock2, group.by = "layer_guess", assay = "Assay4", spatial = "main_reg", alpha = 0.5)

####
## visualize all ####
####

layer_guess <- c("Layer3", "Layer1", "WM", "Layer5", "Layer6", "Layer2", "Layer4")
colors <- setNames(object = hue_pal(7), nm = layer_guess)
colors <- as.list(colors)
g1 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay1", crop = TRUE, colors = colors) |>
  addSpatialLayer(SRBlock, group.by = "layer_guess", assay = "Assay2", alpha = 0.5, spatial = "main", colors = colors) + 
  theme(legend.position = "none") + labs(title = "")
g2 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay1", crop = TRUE, spatial = "main_reg", colors = colors) |>
  addSpatialLayer(SRBlock, group.by = "layer_guess", assay = "Assay2", alpha = 0.5, spatial = "main_reg", colors = colors) +
  theme(legend.position = "none") + labs(title = "") + scale_y_reverse()
g3 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay3", crop = TRUE, colors = colors) |>
  addSpatialLayer(SRBlock, group.by = "layer_guess", assay = "Assay4", alpha = 0.5, spatial = "main", colors = colors) +
  theme(legend.position = "none") + labs(title = "")
g4 <- vrSpatialPlot(SRBlock, group.by = "layer_guess", assay = "Assay3", crop = TRUE, spatial = "main_reg", colors = colors) |>
  addSpatialLayer(SRBlock, group.by = "layer_guess", assay = "Assay4", alpha = 0.5, spatial = "main_reg", colors = colors) +
  theme(legend.position = "none") + labs(title = "") + scale_y_reverse()
g5 <- vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay1", crop = TRUE, colors = colors) |>
  addSpatialLayer(SRBlock2, group.by = "layer_guess", assay = "Assay2", alpha = 0.5, spatial = "main", colors = colors) + 
  theme(legend.position = "none") + labs(title = "")
g6 <- vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay1", crop = TRUE, spatial = "main_reg", colors = colors) |>
  addSpatialLayer(SRBlock2, group.by = "layer_guess", assay = "Assay2", alpha = 0.5, spatial = "main_reg", colors = colors) + 
  theme(legend.position = "none") + labs(title = "") + scale_y_reverse()
g7 <- vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay3", crop = TRUE, colors = colors) |>
  addSpatialLayer(SRBlock2, group.by = "layer_guess", assay = "Assay4", alpha = 0.5, spatial = "main", colors = colors) +
  theme(legend.position = "none") + labs(title = "")
g8 <- vrSpatialPlot(SRBlock2, group.by = "layer_guess", assay = "Assay3", crop = TRUE, spatial = "main_reg", colors = colors) |>
  addSpatialLayer(SRBlock2, group.by = "layer_guess", assay = "Assay4", alpha = 0.5, spatial = "main_reg", colors = colors) +
  theme(legend.position = "none") + labs(title = "") + scale_y_reverse()
colors <- colors[unique(SRBlock3$layer_guess)]
g9 <- vrSpatialPlot(SRBlock3, group.by = "layer_guess", assay = "Assay1", crop = TRUE, colors = colors) |>
  addSpatialLayer(SRBlock3, group.by = "layer_guess", assay = "Assay2", alpha = 0.5, spatial = "main", colors = colors) +
  theme(legend.position = "none") + labs(title = "")
g10 <- vrSpatialPlot(SRBlock3, group.by = "layer_guess", assay = "Assay1", crop = TRUE, spatial = "main_reg", colors = colors) |>
  addSpatialLayer(SRBlock3, group.by = "layer_guess", assay = "Assay2", alpha = 0.5, spatial = "main_reg", colors = colors) + 
  theme(legend.position = "none") + labs(title = "") + scale_y_reverse()
g11 <- vrSpatialPlot(SRBlock3, group.by = "layer_guess", assay = "Assay3", crop = TRUE, colors = colors) |>
  addSpatialLayer(SRBlock3, group.by = "layer_guess", assay = "Assay4", alpha = 0.5, spatial = "main", colors = colors) +
  theme(legend.position = "none") + labs(title = "")
g12 <- vrSpatialPlot(SRBlock3, group.by = "layer_guess", assay = "Assay3", crop = TRUE, spatial = "main_reg", colors = colors) |>
  addSpatialLayer(SRBlock3, group.by = "layer_guess", assay = "Assay4", alpha = 0.5, spatial = "main_reg", colors = colors) + 
  theme(legend.position = "none") + labs(title = "") + scale_y_reverse()

# ggpubr::ggarrange(plotlist = list(g1,g3,g5,g7,g9,g11,g2,g4,g6,g8,g10,g12), nrow = 2, ncol = 6)
ggpubr::ggarrange(plotlist = list(g2,g4,g6,g8,g10,g12), nrow = 2, ncol = 6)
ggsave("results/paste_alloverlay.pdf", width = 36, height = 7)