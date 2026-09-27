# library(VoltRon)
library(ggpubr)
library(magick)
library(patchwork)
library(ggplot2)
library(sesplot) # remotes::install_github("smin95/sesplot")
####
# Xenium vs Visium ####
####

# library(VoltRon)
Xen_R1 <- importXenium("../../../../../data/10X_Xenium_Visium/Xenium_R1/outs/", sample_name = "XeniumR1",
                       resolution_level = 7, overwrite_resolution = TRUE, import_molecules = FALSE)
# Xen_R2 <- importXenium("../../../../../data/10X_Xenium_Visium/Xenium_R2/outs/", sample_name = "XeniumR2",
#                        resolution_level = 3, overwrite_resolution = TRUE, import_molecules = FALSE)
Vis <- importVisium("../../../../../data/10X_Xenium_Visium/Visium/", sample_name = "VisiumR1", resolution_level = "hires")

# insert STAlign coordinates
new_coords <- read.csv("results/Visium_vs_XeniumR1/xenium_to_visium_wsireg.csv.gz")
new_coords <- new_coords[,c("x", "y")]
new_coords <- new_coords * (0.09271707 / 0.02781512)
rownames(new_coords) <- vrSpatialPoints(Xen_R1)
imageinfo <- magick::image_info(vrImages(Vis))
new_coords$y <- imageinfo$width - new_coords$y
new_coords$y <- new_coords$y - imageinfo$width/2
new_coords$x <- new_coords$x - imageinfo$height/2
R <- matrix(c(0,-1,1,0), ncol = 2, byrow = TRUE)
new_coords_rot <- as.matrix(new_coords) %*% t(R)
colnames(new_coords_rot) <- colnames(new_coords)
new_coords_rot <- as.data.frame(new_coords_rot)
new_coords_rot$y <- new_coords_rot$y + imageinfo$height/2
new_coords_rot$x <- new_coords_rot$x + imageinfo$width/2
# new_coords_rot <- new_coords
# new_coords_rot$y <- imageinfo$width - new_coords_rot$y
vrImages(Xen_R1[["Assay1"]], reg = TRUE) <- formImage(coords = new_coords_rot)
vrMainSpatial(Xen_R1[["Assay1"]]) <- "main_reg"
# new_coords <- read.csv("results/Visium_vs_XeniumR1/xenium_to_visium_wsireg_nl.csv.gz")
# new_coords <- new_coords[,c("x", "y")]
# rownames(new_coords) <- vrSpatialPoints(Xen_R1)
# imageinfo <- magick::image_info(vrImages(Vis))
# new_coords$y <- imageinfo$height - new_coords$y
# vrImages(Xen_R1[["Assay1"]], name = "main_reg2") <- formImage(coords = new_coords)
# vrMainSpatial(Xen_R1[["Assay1"]]) <- "main_reg2"

merge_list <- list(Xen_R1, Vis)
VRBlock <- merge(merge_list[[1]], merge_list[-1], samples = "DLPFC_Block")
VRBlock

####
# visualize ####
####

set.seed(1)
td1 <- tempfile(fileext = ".png")
g1 <- vrSpatialPlot(VRBlock, assay = "Assay2", alpha = 0)
ggsave(filename = td1, plot = g1, width = 8, height = 8)
td2 <- tempfile(fileext = ".png")
g2 <- vrSpatialPlot(VRBlock, assay = "Assay2", alpha = 0) |>
  addSpatialLayer(VRBlock, assay = "Assay1", spatial = "main_reg")
ggsave(filename = td2, plot = g2, width = 8, height = 8)
c(image_read(td1), image_read(td2)) %>%
  magick::image_join()
# 
# set.seed(1)
# td1 <- tempfile(fileext = ".png")
# g1 <- vrSpatialPlot(VRBlock, assay = "Assay2", alpha = 0)
# ggsave(filename = td1, plot = g1, width = 8, height = 8)
# td2 <- tempfile(fileext = ".png")
# g2 <- vrSpatialPlot(VRBlock, assay = "Assay2", alpha = 0) |> 
#   addSpatialLayer(VRBlock, assay = "Assay3", spatial = "main_reg")
# ggsave(filename = td2, plot = g2, width = 8, height = 8)
# c(image_read(td1), image_read(td2)) %>% 
#   magick::image_join()

####
# compare spots and pseudospots ####
####

vrMainSpatial(Xen_R1[["Assay1"]]) <- "main_reg"
VRBlock <- transferData(VRBlock, from = "Assay1", to = "Assay2", new_feature_name = "RNA_pseudoXenium1")
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
# p2 <- ggplot(data = datax, mapping = aes(x = Visium, y = Xenium)) + 
#   geom_point(color = "darkgreen") + 
#   theme_classic() + labs(title = "Specificity (TACSTD2)")
ggplot(data = datax, mapping = aes(x = Visium, y = Xenium)) +
  geom_point(shape = 21, fill = "grey60", color = "grey60", size = 3) +
  sm_statCorr(
    color = "black", corr_method = "pearson", R2 = TRUE
  ) + 
  theme_classic() + labs(title = "Specificity (TACSTD2)")
ggsave(filename = "../../../../../Documents/VoltRon Paper/Nature Methods Revision/Images/Supplementary Material/ImageAlign/images/wsireg_scatterplot_TACSTD2.pdf", plot = last_plot(), 
       width = 7, height = 7, device = "pdf")

# calculate lee's l measure
vrMainAssay(VRBlock) <- "Visium"
VRBlock <- getSpatialNeighbors(VRBlock, method = "radius")
matgraph <- igraph::as_adjacency_matrix(vrGraph(VRBlock))
matgraph <- matgraph[names(data2),names(data2)]
neigh <- spdep::mat2listw(matgraph)
lee.test(x = datax[,1], 
         y = datax[,2], 
         listw = neigh, zero.policy=attr(neigh, "zero.policy"),
         alternative="greater")
# Lee's L statistic standard deviate = 5.0777, p-value = 1.91e-07
# alternative hypothesis: greater
# sample estimates:
# Lee's L statistic       Expectation          Variance 
# 2.477498e-02      2.580604e-03      1.910538e-05 

### ERBB2 ####

# correlation plot of ERBB2 across two assays
data1 <- vrData(VRBlock, assay = "Assay2", features = "ERBB2", feat_type = "RNA")
data2 <- vrData(VRBlock, assay = "Assay2", features = "ERBB2", feat_type = "RNA_pseudoXenium1")
data1 <- data1[1,,drop=TRUE]
data2 <- data2[1,,drop=TRUE]
datax <- data.frame(Visium = log10(data1[names(data2)] + 1), 
                    Xenium = log10(data2 + 1))
# p2 <- ggplot(data = datax, mapping = aes(x = Visium, y = Xenium)) + 
#   geom_point(color = "darkgreen") + 
#   theme_classic() + labs(title = "Specificity (ERBB2)")
ggplot(data = datax, mapping = aes(x = Visium, y = Xenium)) +
  geom_point(shape = 21, fill = "grey60", color = "grey60", size = 3) +
  sm_statCorr(
    color = "black", corr_method = "pearson", R2 = TRUE
  ) + 
  theme_classic() + labs(title = "Specificity (ERBB2)")
ggsave(filename = "../../../../../Documents/VoltRon Paper/Nature Methods Revision/Images/Supplementary Material/ImageAlign/images/wsireg_scatterplot_ERBB2.pdf", plot = last_plot(), 
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
# Lee's L statistic standard deviate = 35.529, p-value < 2.2e-16
# alternative hypothesis: greater
# sample estimates:
# Lee's L statistic       Expectation          Variance 
# 1.855978e-01      3.021969e-02      1.912503e-05 

## Visium vs XeniumR2 ####

# # correlation plot of all genes across two assays
# data1 <- vrData(VRBlock, assay = "Assay2", feat_type = "RNA")
# data2 <- vrData(VRBlock, assay = "Assay2", feat_type = "RNA_pseudoXenium2")
# data1 <- rowMeans(data1)
# data2 <- rowMeans(data2)
# datax <- data.frame(Visium = log10(data1[names(data2)]),
#                     Xenium = log10(data2))
# p3 <- ggplot(data = datax, mapping = aes(x = Visium, y = Xenium)) + 
#   geom_point(color = "darkgreen") + 
#   geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") + 
#   geom_abline(slope = 1, intercept = 1, color = "darkgreen") +
#   theme_classic() + labs(title = "Sensitivity (all panel genes)")
# 
# # correlation plot of TACSTD2 across two assays
# data1 <- vrData(VRBlock, assay = "Assay2", features = "TACSTD2", feat_type = "RNA")
# data2 <- vrData(VRBlock, assay = "Assay2", features = "TACSTD2", feat_type = "RNA_pseudoXenium2")
# data1 <- data1[1,,drop=TRUE]
# data2 <- data2[1,,drop=TRUE]
# datax <- data.frame(Visium = log10(data1[names(data2)] + 1), 
#                     Xenium = log10(data2 + 1))
# p4 <- ggplot(data = datax, mapping = aes(x = Visium, y = Xenium)) + 
#   geom_point(color = "darkgreen") + 
#   theme_classic() + labs(title = "Specificity (TACSTD2)")

## visualize all ####
# (p1 | p3) / (p2 | p4)
# ggsave(filename = "../../allresults/BreastCancer/figures/scatterplot/stalign_XenvsVisium.pdf", 
#        plot = p1 | p2 | p3 | p4,
#        device = "pdf", width = 24, height = 6)
(p1) | (p2)
ggsave(filename = "../../allresults/BreastCancer/figures/scatterplot/wsireg_XenvsVisium.pdf", 
       plot = p1 | p2,
       device = "pdf", width = 12, height = 6)

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
ggsave(filename = "../../allresults/BreastCancer/figures/visium_align_zoom/wsireg_subset_breastcancervisium.pdf", 
       plot = last_plot(), width = 7, height = 7)

## visualize plot overlay ####
img1 <- vrImages(VRBlock, assay = "Assay2")
points <- vrCoordinates(VRBlock, assay = "Assay1")
points <- as.data.frame(points[,c(1,2)])
# points <- points * (0.17011142 / 0.051033426)
colnames(points) <- c("x", "y")
# points$y <- magick::image_info(img1)$height - points$y
# img1 <- image_resize(img1, geometry = magick::geometry_size_percent(20))
# img1 <- image_rotate(img1, degrees = 90)
# points <- points*0.20
g1 <- magick::image_ggplot(img1) + 
  stat_bin_2d(data = points, mapping = aes(x = x, y = y), color = "black", size = 0.05, alpha = 0.5, bins = 200)
datax <- g1@layers$stat_bin_2d$data
magick::image_ggplot(img1) + 
  geom_point(data = datax, mapping = aes(x = x, y = y), color = "black", size = 0.001, alpha = 0.5)
ggsave(filename = "../../allresults/BreastCancer/figures/visium_overlay/wsireg_XenvsVisium_overlay.png", 
       device = "png", plot = last_plot(), width = 6, height = 8)

####
# Lee's test across all genes ####
####

VRBlock <- transferData(VRBlock, from = "Assay1", to = "Assay2", new_feature_name = "RNA_pseudoXenium1")
vrMainFeatureType(VRBlock[["Assay2"]]) <- "RNA"
features1 <- vrFeatures(VRBlock)
vrMainFeatureType(VRBlock[["Assay2"]]) <- "RNA_pseudoXenium1"
features2 <- vrFeatures(VRBlock)
selected_features <- intersect(features1, features2)

# Lee's Test
leesstats <- NULL
for(feat in selected_features){
  
  # correlation plot of TACSTD2 across two assays
  data1 <- vrData(VRBlock, assay = "Assay2", features = feat, feat_type = "RNA")
  data2 <- vrData(VRBlock, assay = "Assay2", features = feat, feat_type = "RNA_pseudoXenium1")
  data1 <- data1[1,,drop=TRUE]
  data2 <- data2[1,,drop=TRUE]
  datax <- data.frame(Visium = log10(data1[names(data2)] + 1), 
                      Xenium = log10(data2 + 1))
  vrMainAssay(VRBlock) <- "Visium"
  VRBlock <- getSpatialNeighbors(VRBlock, method = "radius")
  matgraph <- igraph::as_adjacency_matrix(vrGraph(VRBlock))
  matgraph <- matgraph[names(data2),names(data2)]
  neigh <- spdep::mat2listw(matgraph, )
  res <- lee.test(x = datax[,1], 
                  y = datax[,2], 
                  listw = neigh, zero.policy=attr(neigh, "zero.policy"),
                  alternative="greater") 
  
  # add stats
  res <- res$estimate
  res <- c("gene" = feat, res)
  print(res)
  leesstats <- rbind(leesstats, res)
}
