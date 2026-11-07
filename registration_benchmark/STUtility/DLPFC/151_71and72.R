# library(VoltRon)
library(patchwork)
library(anndataR)
library(anndata)
library(STutility)
library(patchwork)
library(ggpubr)
# source("../../allresults/DLFPC/accuracy.R")

set.seed(1)

####
# import ####
####

DLPFC_1_dir <- "../../../../../data/10X_Visium_DLPFC/151671//"
DLPFC_2_dir <- "../../../../../data/10X_Visium_DLPFC/151672/"
dir_list <- c(DLPFC_1_dir, DLPFC_2_dir)
dir_list <- sapply(dir_list, function(x) {
  tmp <- c("filtered_feature_bc_matrix.h5", "spatial/tissue_positions_list.txt", "spatial/tissue_hires_image.png", "spatial/scalefactors_json.json")
  tmp <- sapply(tmp, function(y) file.path(x, y))
  print(tmp)
  names(tmp) <- c("samples", "spotfiles", "imgs", "json")
  tmp
})
dir_list <- as.data.frame(t(dir_list))
se <- InputFromTable(infotable = dir_list, 
                     platform =  "Visium")

# layers
layer_guess1 <- readRDS("../../../benchmarks/allresults/DLFPC/151671_labels.rds")
layer_guess2 <- readRDS("../../../benchmarks/allresults/DLFPC/151672_labels.rds")
names(layer_guess1) <- paste0(names(layer_guess1), "_1")
names(layer_guess2) <- paste0(names(layer_guess2), "_2")
layer_guess <- c(layer_guess1, layer_guess2)
layer_guess <- setNames(as.character(layer_guess), names(layer_guess))
se$layer_guess <- layer_guess[Cells(se)]
tmp <- se$layer_guess
tmp[is.na(tmp)] <- "none"
se$layer_guess <- tmp

# Process
se <- NormalizeData(se, scale.factor = 10000)

# visualize
layout(matrix(1:3,nrow =3))
se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)
ImagePlot(se, method = "raster", type = "raw", ncols = 2)

# mask
se2 <- MaskImages(object = se)
ImagePlot(se2, ncols = 2, method = "raster", type = "masked")
ImagePlot(se2, ncols = 2, method = "raster", type = "masked.masks", darken = TRUE)

# align
se3 <- AlignImages(se2)
ImagePlot(se3, method = "raster", type = "processed", ncols = 2)
ImagePlot(se3, method = "raster", type = "processed.masks", ncols = 2)

# visualize
# ST.FeaturePlot(se3, features = "layer_guess", ncol = 2, indices = 1)
# ST.FeaturePlot(se, features = "orig.ident", ncol = 2)

####
# visualize ####
####

####
## 1 vs 2 ####
####

set.seed(1)
td1 <- tempfile(fileext = ".png")
colors <- setNames(object = hue_pal(8), nm = unique(se$layer_guess))
metadata <- cbind(se3@tools$Staffli@meta.data, se3@meta.data)
metadata <- metadata[,c("warped_x", "warped_y", "sample", "layer_guess")]
colnames(metadata) <- c("x", "y", "sample", "layer_guess")
metadata1 <- metadata[metadata$sample==1,c("x", "y", "layer_guess")]
limits <- apply(metadata1[,1:2], 2, range) + c(-50,50)
limits
g1 <- ggplot(metadata1, aes(x=x,y=y, color=layer_guess)) + 
  geom_point() + 
  theme_pubclean() + NoLegend() + 
  xlim(limits[1,1], limits[2,1]) + ylim(limits[1,2], limits[2,2])
ggsave(filename = td1, plot = g1, width = 8, height = 8)
td2 <- tempfile(fileext = ".png")
metadata2 <- metadata[metadata$sample==2,c("x", "y", "layer_guess")]
g2 <- ggplot(metadata2, aes(x=x,y=y, color=layer_guess)) + 
  geom_point() + 
  theme_pubclean() + NoLegend() + 
  xlim(limits[1,1], limits[2,1]) + ylim(limits[1,2], limits[2,2])
ggsave(filename = td2, plot = g2, width = 8, height = 8)
c(image_read(td1), image_read(td2)) %>% 
  magick::image_join()
write.table(metadata2, "results/151672_aligned.rds")

