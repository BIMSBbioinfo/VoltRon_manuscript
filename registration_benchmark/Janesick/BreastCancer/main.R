# library(VoltRon)
library(ggpubr)
library(magick)
library(patchwork)
library(jsonlite)
library(ggplot2)

####
# Xenium vs Visium ####
####

# library(VoltRon)
Xen_R1 <- importXenium("../../../../../data/10X_Xenium_Visium/Xenium_R1/outs/", sample_name = "XeniumR1",
                       resolution_level = 7, overwrite_resolution = TRUE, import_molecules = FALSE)
Vis <- importVisium("../../../../../data/10X_Xenium_Visium/Visium/", sample_name = "VisiumR1", resolution_level = "hires")

# insert Janesick transformation matrix
new_coords <- read.csv("results/xenium_cell_centroids_visium_high_res.csv", row.names = 1)
colnames(new_coords) <- c("x", "y")
rownames(new_coords) <- vrSpatialPoints(Xen_R1)
imageinfo <- magick::image_info(vrImages(Vis))
new_coords$y <- imageinfo$height - new_coords$y
vrImages(Xen_R1[["Assay1"]], reg = TRUE) <- formImage(coords = new_coords)
vrMainSpatial(Xen_R1[["Assay1"]]) <- "main_reg"
merge_list <- list(Xen_R1, Vis)
VRBlock <- merge(merge_list[[1]], merge_list[-1], samples = "DLPFC_Block")
VRBlock

####
# visualize ####
####

####
## 1 vs 2 ####
####

# set.seed(1)
# td1 <- tempfile(fileext = ".png")
# g1 <- vrSpatialPlot(VRBlock, assay = "Assay2", alpha = 0)
# ggsave(filename = td1, plot = g1, width = 8, height = 8)
# td2 <- tempfile(fileext = ".png")
# g2 <- vrSpatialPlot(VRBlock, assay = "Assay2", alpha = 0) |> 
#   addSpatialLayer(VRBlock, assay = "Assay1", spatial = "main_reg")
# ggsave(filename = td2, plot = g2, width = 8, height = 8)
# c(image_read(td1), image_read(td2)) %>% 
#   magick::image_join()

####
# compare spots and pseudospots ####
####

VRBlock <- transferData(VRBlock, from = "Assay1", to = "Assay2", new_feature_name = "RNA_pseudoXenium1")
# VRBlock <- transferData(VRBlock, from = "Assay3", to = "Assay2", new_feature_name = "RNA_pseudoXenium2")
vrMainFeatureType(VRBlock[["Assay2"]]) <- "RNA_pseudoXenium1"
vrMainFeatureType(VRBlock, assay = "all")

## Visium vs XeniumR1 ####

# correlation plot of all genes across two assays
data1 <- vrData(VRBlock, assay = "Assay2", feat_type = "RNA")
data2 <- vrData(VRBlock, assay = "Assay2", feat_type = "RNA_pseudoXenium1")
data1 <- rowMeans(data1)
data2 <- rowMeans(data2)
datax <- data.frame(Visium = log10(data1[names(data2)]),
                    Xenium = log10(data2))
p1 <- ggplot(data = datax, mapping = aes(x = Visium, y = Xenium)) + 
  geom_point(color = "darkgreen") + 
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") + 
  geom_abline(slope = 1, intercept = 1, color = "darkgreen") +
  theme_classic() + labs(title = "Sensitivity (all panel genes)")

### TACSTD2 ####

# correlation plot of TACSTD2 across two assays
data1 <- vrData(VRBlock, assay = "Assay2", features = "TACSTD2", feat_type = "RNA")
data2 <- vrData(VRBlock, assay = "Assay2", features = "TACSTD2", feat_type = "RNA_pseudoXenium1")
data1 <- data1[1,,drop=TRUE]
data2 <- data2[1,,drop=TRUE]
datax <- data.frame(Visium = log10(data1[names(data2)] + 1), 
                    Xenium = log10(data2 + 1))
ggplot(data = datax, mapping = aes(x = Visium, y = Xenium)) +
  geom_point(shape = 21, fill = "grey60", color = "grey60", size = 3) +
  sm_statCorr(
    color = "black", corr_method = "pearson", R2 = TRUE
  ) + 
  theme_classic() + labs(title = "Specificity (TACSTD2)")
ggsave(filename = "../../../../../Documents/VoltRon Paper/Nature Methods Revision/Images/Supplementary Material/ImageAlign/images/janesick_scatterplot_TACSTD2.pdf", plot = last_plot(), 
       width = 7, height = 7, device = "pdf")

# calculate lee's l measure
vrMainAssay(VRBlock) <- "Visium"
VRBlock <- getSpatialNeighbors(VRBlock, method = "radius")
matgraph <- igraph::as_adjacency_matrix(vrGraph(VRBlock))
matgraph <- matgraph[names(data2),names(data2)]
neigh <- spdep::mat2listw(matgraph, )
lee.test(x = datax[,1], 
         y = datax[,2], 
         listw = neigh, zero.policy=attr(neigh, "zero.policy"),
         alternative="greater")
# Lee's L statistic standard deviate = 110.78, p-value < 2.2e-16
# alternative hypothesis: greater
# sample estimates:
# Lee's L statistic       Expectation          Variance 
# 6.620825e-01      1.443742e-01      2.183901e-05 

### ERBB2 ####

# correlation plot of ERBB2 across two assays
data1 <- vrData(VRBlock, assay = "Assay2", features = "ERBB2", feat_type = "RNA")
data2 <- vrData(VRBlock, assay = "Assay2", features = "ERBB2", feat_type = "RNA_pseudoXenium1")
data1 <- data1[1,,drop=TRUE]
data2 <- data2[1,,drop=TRUE]
datax <- data.frame(Visium = log10(data1[names(data2)] + 1), 
                    Xenium = log10(data2 + 1))
ggplot(data = datax, mapping = aes(x = Visium, y = Xenium)) +
  geom_point(shape = 21, fill = "grey60", color = "grey60", size = 3) +
  sm_statCorr(
    color = "black", corr_method = "pearson", R2 = TRUE
  ) + 
  theme_classic() + labs(title = "Specificity (ERBB2)")
ggsave(filename = "../../../../../Documents/VoltRon Paper/Nature Methods Revision/Images/Supplementary Material/ImageAlign/images/janesick_scatterplot_ERBB2.pdf", plot = last_plot(), 
       width = 7, height = 7, device = "pdf")

# calculate lee's l measure
vrMainAssay(VRBlock) <- "Visium"
VRBlock <- getSpatialNeighbors(VRBlock, method = "radius")
matgraph <- igraph::as_adjacency_matrix(vrGraph(VRBlock))
matgraph <- matgraph[names(data2),names(data2)]
neigh <- spdep::mat2listw(matgraph, )
lee.test(x = datax[,1], 
         y = datax[,2], 
         listw = neigh, zero.policy=attr(neigh, "zero.policy"),
         alternative="greater")

####
# compare visualization ####
####

VRBlock_assay1 <- subset(VRBlock, assay = "Assay1")
vrImages(VRBlock_assay1[["Assay1"]], name = "main_reg", channel = "H&E") <- vrImages(VRBlock, assay = "Assay2")
vrMainChannel(VRBlock_assay1[["Assay1"]]) <- "H&E"
VRBlock_assay1 <- subset(VRBlock_assay1, image = readRDS("../../allresults/BreastCancer/figures/visium_align_zoom/Visium_vis_subset.rds"))
vrSpatialPlot(VRBlock_assay1, spatial = "main_reg", channel = "H&E", color = list(DLPFC_Block = "black")) + 
  theme(legend.position = "none") + 
  labs(title = "")
ggsave(filename = "../../allresults/BreastCancer/figures/visium_align_zoom/Janesick_subset_breastcancervisium.pdf", 
       plot = last_plot(), width = 7, height = 7)

## visualize plot overlay ####
img1 <- vrImages(VRBlock, assay = "Assay2")
points <- vrCoordinates(VRBlock, assay = "Assay1")
points <- as.data.frame(points[,c(1,2)])
colnames(points) <- c("x", "y")
# points$y <- magick::image_info(img1)$height - points$y
img1 <- image_resize(img1, geometry = magick::geometry_size_percent(20))
points <- points*0.20
g1 <- magick::image_ggplot(img1) + 
  stat_bin_2d(data = points, mapping = aes(x = x, y = y), color = "black", size = 0.05, alpha = 0.5, bins = 200)
datax <- g1@layers$stat_bin_2d$data
magick::image_ggplot(img1) + 
  geom_point(data = datax, mapping = aes(x = x, y = y), color = "black", size = 0.05, alpha = 0.5)
ggsave(filename = "../../allresults/BreastCancer/figures/visium_overlay/janesick_XenvsVisium_overlay.pdf", 
       device = "png", plot = last_plot(), width = 6, height = 8)