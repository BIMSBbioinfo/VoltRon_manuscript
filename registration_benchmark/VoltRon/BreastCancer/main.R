# library(VoltRon)
library(ggpubr)
library(magick)
library(patchwork)
library(ggplot2)

####
# tests ####
####

# import data
XeniumR1_image <- magick::image_read(paste0("../../../../ImageAlignmentBenchmark/data/BreastCancer/Xenium_R1/Xenium/image_3.tiff"))
DAPI_image <- magick::image_read(paste0("../../../../ImageAlignmentBenchmark/data/BreastCancer/Xenium_R1/IF/DAPI.tiff"))
HE_image <- magick::image_read(paste0("../../../../ImageAlignmentBenchmark/data/BreastCancer/Xenium_R1/HE/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image_highres.tif"))

####
# XeniumR1 vs HE ####
####

####
## get parameters ####
####

ref_image <- HE_image
query_image <- XeniumR1_image
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(query_image), 
                               mapping_parameters = readRDS("parameters/XeniumR1_vs_HE_homography.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/XeniumR1_vs_HE_homography.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(query_image), 
                               mapping_parameters = readRDS("parameters/XeniumR1_vs_HE_homography_TPS.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/XeniumR1_vs_HE_homography_TPS.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(query_image), 
                               mapping_parameters = readRDS("parameters/XeniumR1_vs_HE_affine.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/XeniumR1_vs_HE_affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(query_image), 
                               mapping_parameters = readRDS("parameters/XeniumR1_vs_HE_affinetps.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/XeniumR1_vs_HE_affinetps.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(query_image), 
                               mapping_parameters = readRDS("parameters/XeniumR1_vs_HE_manual_Affine.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/XeniumR1_vs_HE_manual_Affine.rds")

####
## transformation of validation landmarks ####
####

# landmarks and images
ref_image <- HE_image
query_image <- magick::image_flip(XeniumR1_image)
HE_landmarks <- read.csv("../../landmarks/BreastCancer/XeniumR1_vs_HE/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image_highres_landmarks.csv", row.names = 1) 
Xenium_landmarks <- read.csv("../../landmarks/BreastCancer/XeniumR1_vs_HE/Xenium_R1_image_3.csv", row.names = 1) 
imageinfo_query <- magick::image_info(query_image)
imageinfo_ref <- magick::image_info(ref_image)

# transform validation markers (auto)
mapping_parameters <- readRDS("parameters/XeniumR1_vs_HE_homography.rds")
Xenium_landmarks_auto <- Xenium_landmarks[,c("x","y")]
# Xenium_landmarks_auto$y <- imageinfo_query$height - Xenium_landmarks_auto$y
Xenium_landmarks_auto[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_auto), mapping_parameters$mapping$`2`)
Xenium_landmarks_auto$y <- imageinfo_ref$height - Xenium_landmarks_auto$y
diff_auto <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_auto[,c("x", "y")])^2))

# transform validation markers (auto + TPS)
mapping_parameters <- readRDS("parameters/XeniumR1_vs_HE_homography_TPS.rds")
Xenium_landmarks_autotps <- Xenium_landmarks[,c("x","y")]
# Xenium_landmarks_autotps$y <- imageinfo_query$height - Xenium_landmarks_autotps$y
Xenium_landmarks_autotps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_autotps), mapping_parameters$mapping$`2`)
Xenium_landmarks_autotps$y <- imageinfo_ref$height - Xenium_landmarks_autotps$y
diff_autotps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_autotps[,c("x", "y")])^2))

# transform validation markers (affine)
mapping_parameters <- readRDS("parameters/XeniumR1_vs_HE_affine.rds")
Xenium_landmarks_affine <- Xenium_landmarks[,c("x","y")]
# Xenium_landmarks_affine$y <- imageinfo_query$height - Xenium_landmarks_affine$y
Xenium_landmarks_affine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine$y <- imageinfo_ref$height - Xenium_landmarks_affine$y
diff_affine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine[,c("x", "y")])^2))

# transform validation markers (affine + TPS)
mapping_parameters <- readRDS("parameters/XeniumR1_vs_HE_homography.rds")
Xenium_landmarks_affinetps <- Xenium_landmarks[,c("x","y")]
# Xenium_landmarks_affinetps$y <- imageinfo_query$height - Xenium_landmarks_affinetps$y
Xenium_landmarks_affinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_affinetps$y <- imageinfo_ref$height - Xenium_landmarks_affinetps$y
diff_affinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affinetps[,c("x", "y")])^2))

# transform validation markers (Manual affine)
mapping_parameters <- readRDS("parameters/XeniumR1_vs_HE_manual_Affine.rds")
Xenium_landmarks_manualaffine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffine$y <- imageinfo_query$height - Xenium_landmarks_manualaffine$y
Xenium_landmarks_manualaffine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffine), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffine$y <- imageinfo_ref$height - Xenium_landmarks_manualaffine$y
diff_manualaffine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffine[,c("x", "y")])^2))

# transform validation markers (affine + BSpline)
mapping <- readRDS("parameters/XeniumR1_vs_HE_affine.rds")
mapping$Method <- "Affine + Non-Rigid"
mapping$nonrigid <- "BSpline (SimpleITK)"
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(XeniumR1_image),
                               mapping_parameters = mapping, interactive = TRUE)
mapping_parameters <- xen_reg$mapping_parameters
Xenium_landmarks_affine_nr <- Xenium_landmarks[,c("x","y")]
# Xenium_landmarks_affine_nr$y <- imageinfo_query$height - Xenium_landmarks_affine_nr$y
Xenium_landmarks_affine_nr[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine_nr), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine_nr$y <- imageinfo_ref$height - Xenium_landmarks_affine_nr$y
diff_affine_nr <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine_nr[,c("x", "y")])^2))

# compare markers
# results <- data.frame(HE_landmarks[,c("x","y")], diff_manualtps, diff_manualhomotps, diff_auto, diff_autotps)
# colnames(results) <- c('x_landmark', 'y_landmark', "diff_manualtps", "diff_manualhomotps", "diff_auto", "diff_autotps")
results <- data.frame(HE_landmarks[,c("x","y")], diff_affine, diff_affinetps, diff_auto, diff_autotps, diff_manualaffine, diff_affine_nr)
colnames(results) <- c('x_landmark', 'y_landmark', "diff_affine", "diff_affinetps", "diff_auto", "diff_autotps", "diff_manualaffine", "diff_affine_nr")
dir.create("results/XeniumR1_vs_HE/")
write.csv(results, file = "results/XeniumR1_vs_HE/all_landmarks.csv", row.names = FALSE, quote = FALSE)

####
## visualize and check ####
####

# visualize
g1 <- image_ggplot(image_resize(XeniumR1_image, geometry = magick::geometry_size_percent(20))) +
  geom_point(mapping = aes(x = x, y = y), data = Xenium_landmarks[,c("x","y")]/5, colour = "blue", size = 3)

# visualize
g2 <- image_ggplot(image_resize(HE_image, geometry = magick::geometry_size_percent(20))) +
  geom_point(mapping = aes(x = x, y = y), data = HE_landmarks[,c("x","y")]/5, colour = "blue", size = 3) + 
  geom_point(mapping = aes(x = x, y = y), data = Xenium_landmarks_auto[,c("x","y")]/5, colour = "red", size = 3)

####
# XeniumR1 vs DAPI ####
####

####
## get parameters ####
####

ref_image <- DAPI_image
query_image <- XeniumR1_image
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(query_image), 
                               mapping_parameters = readRDS("parameters/XeniumR1_vs_DAPI_homography.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/XeniumR1_vs_DAPI_homography.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(query_image), 
                               mapping_parameters = readRDS("parameters/XeniumR1_vs_DAPI_homography.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/XeniumR1_vs_HE_homography_TPS.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(query_image), 
                               mapping_parameters = readRDS("parameters/XeniumR1_vs_DAPI_affine.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/XeniumR1_vs_DAPI_affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(query_image), 
                               mapping_parameters = readRDS("parameters/XeniumR1_vs_DAPI_affinetps.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/XeniumR1_vs_DAPI_affinetps.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(query_image), 
                               mapping_parameters = readRDS("parameters/XeniumR1_vs_DAPI_manual_Affine.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/XeniumR1_vs_DAPI_manual_Affine.rds")

####
## transformation of validation landmarks ####
####

# landmarks and images
ref_image <- DAPI_image
query_image <- magick::image_rotate(XeniumR1_image, degrees = 180)
DAPI_landmarks <- read.csv("../../landmarks/BreastCancer/XeniumR1_vs_DAPI/Xenium_DAPI_landmarks.csv", row.names = 1) 
Xenium_landmarks <- read.csv("../../landmarks/BreastCancer/XeniumR1_vs_DAPI/Xenium_R1_image_3.csv", row.names = 1) 
imageinfo_query <- magick::image_info(query_image)
imageinfo_ref <- magick::image_info(ref_image)

# transform validation markers (auto)
mapping_parameters <- readRDS("parameters/XeniumR1_vs_DAPI_homography.rds")
Xenium_landmarks_auto <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_auto$x <- imageinfo_query$width - Xenium_landmarks_auto$x
Xenium_landmarks_auto[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_auto), mapping_parameters$mapping$`2`)
Xenium_landmarks_auto$y <- imageinfo_ref$height - Xenium_landmarks_auto$y
diff_auto <- sqrt(rowSums((DAPI_landmarks[,c("x","y")] - Xenium_landmarks_auto[,c("x", "y")])^2))

# transform validation markers (auto + TPS)
# mapping_parameters <- readRDS("parameters/XeniumR1_vs_DAPI_homography_TPS.rds")
# Xenium_landmarks_autotps <- Xenium_landmarks[,c("x","y")]
# Xenium_landmarks_autotps$y <- imageinfo_query$height - Xenium_landmarks_autotps$y
# Xenium_landmarks_autotps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_autotps), mapping_parameters$mapping$`2`)
# Xenium_landmarks_autotps$y <- imageinfo_ref$height - Xenium_landmarks_autotps$y
# diff_autotps <- sqrt(rowSums((DAPI_landmarks[,c("x","y")] - Xenium_landmarks_autotps[,c("x", "y")])^2))

# transform validation markers (affine)
mapping_parameters <- readRDS("parameters/XeniumR1_vs_DAPI_affine.rds")
Xenium_landmarks_affine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine$x <- imageinfo_query$width - Xenium_landmarks_affine$x
Xenium_landmarks_affine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine$y <- imageinfo_ref$height - Xenium_landmarks_affine$y
diff_affine <- sqrt(rowSums((DAPI_landmarks[,c("x","y")] - Xenium_landmarks_affine[,c("x", "y")])^2))

# transform validation markers (affine + TPS)
mapping_parameters <- readRDS("parameters/XeniumR1_vs_DAPI_affinetps.rds")
Xenium_landmarks_affinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affinetps$x <- imageinfo_query$width - Xenium_landmarks_affinetps$x
Xenium_landmarks_affinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_affinetps$y <- imageinfo_ref$height - Xenium_landmarks_affinetps$y
diff_affinetps <- sqrt(rowSums((DAPI_landmarks[,c("x","y")] - Xenium_landmarks_affinetps[,c("x", "y")])^2))

# transform validation markers (Manual affine)
mapping_parameters <- readRDS("parameters/XeniumR1_vs_DAPI_manual_Affine.rds")
Xenium_landmarks_manualaffine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffine$y <- imageinfo_query$height - Xenium_landmarks_manualaffine$y
Xenium_landmarks_manualaffine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffine), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffine$y <- imageinfo_ref$height - Xenium_landmarks_manualaffine$y
diff_manualaffine <- sqrt(rowSums((DAPI_landmarks[,c("x","y")] - Xenium_landmarks_manualaffine[,c("x", "y")])^2))

# transform validation markers (affine + BSpline)
mapping <- readRDS("parameters/XeniumR1_vs_DAPI_affine.rds")
mapping$Method <- "Affine + Non-Rigid"
mapping$nonrigid <- "BSpline (SimpleITK)"
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(magick::image_rotate(XeniumR1_image, degrees = 180)),
                               mapping_parameters = mapping, interactive = TRUE)
mapping_parameters <- xen_reg$mapping_parameters
Xenium_landmarks_affine_nr <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine_nr$y <- imageinfo_query$height - Xenium_landmarks_affine_nr$y
Xenium_landmarks_affine_nr[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine_nr), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine_nr$y <- imageinfo_ref$height - Xenium_landmarks_affine_nr$y
diff_affine_nr <- sqrt(rowSums((DAPI_landmarks[,c("x","y")] - Xenium_landmarks_affine_nr[,c("x", "y")])^2))

# compare markers
# results <- data.frame(DAPI_landmarks[,c("x","y")], diff_manualtps, diff_manualhomotps, diff_auto, diff_autotps)
# colnames(results) <- c('x_landmark', 'y_landmark', "diff_manualtps", "diff_manualhomotps", "diff_auto", "diff_autotps")
results <- data.frame(DAPI_landmarks[,c("x","y")], diff_auto, diff_affine, diff_affinetps, diff_manualaffine, diff_affine_nr)
colnames(results) <- c('x_landmark', 'y_landmark', "diff_auto", "diff_affine", "diff_affinetps", "diff_manualaffine", "diff_affine_nr")
dir.create("results/XeniumR1_vs_DAPI/")
write.csv(results, file = "results/XeniumR1_vs_DAPI/all_landmarks.csv", row.names = FALSE, quote = FALSE)

####
## visualize and check ####
####

# visualize
g1 <- image_ggplot(image_resize(XeniumR1_image, geometry = magick::geometry_size_percent(20))) +
  geom_point(mapping = aes(x = x, y = y), data = Xenium_landmarks[,c("x","y")]/5, colour = "blue", size = 3)

# visualize
g2 <- image_ggplot(image_resize(DAPI_image, geometry = magick::geometry_size_percent(20))) +
  geom_point(mapping = aes(x = x, y = y), data = DAPI_landmarks[,c("x","y")]/5, colour = "blue", size = 3)

g2 <- image_ggplot(image_resize(DAPI_image, geometry = magick::geometry_size_percent(20))) +
  geom_point(mapping = aes(x = x, y = y), data = DAPI_landmarks[,c("x","y")]/5, colour = "blue", size = 3) + 
  geom_point(mapping = aes(x = x, y = y), data = Xenium_landmarks_auto[,c("x","y")]/5, colour = "red", size = 3)

####
# Xenium vs Visium ####
####

# library(VoltRon)
Xen_R1 <- importXenium("../../../../../data/10X_Xenium_Visium/Xenium_R1/outs/", sample_name = "XeniumR1",
                       resolution_level = 3, overwrite_resolution = TRUE, import_molecules = FALSE)
Xen_R2 <- importXenium("../../../../../data/10X_Xenium_Visium/Xenium_R2/outs/", sample_name = "XeniumR2",
                       resolution_level = 3, overwrite_resolution = TRUE, import_molecules = FALSE)
Vis <- importVisium("../../../../../data/10X_Xenium_Visium/Visium/", sample_name = "VisiumR1", resolution_level = "hires")

# automated
Xen_R2 <- modulateImage(Xen_R2, brightness = 400)
# xen_reg <- registerSpatialData(object_list = list(Xen_R1, Vis, Xen_R2))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/Xenium_vs_Visium_params.rds")
xen_reg <- registerSpatialData(object_list = list(Xen_R1, Vis, Xen_R2), 
                               mapping_parameters = readRDS("parameters/Xenium_vs_Visium_params.rds"), 
                               interactive = FALSE)
merge_list <- xen_reg$registered_spat
VRBlock <- merge(merge_list[[1]], merge_list[-1], samples = "10XBlock")
VRBlock

## compare spots and pseudospots ####

VRBlock <- transferData(VRBlock, from = "Assay1", to = "Assay2", new_feature_name = "RNA_pseudoXenium1")
VRBlock <- transferData(VRBlock, from = "Assay3", to = "Assay2", new_feature_name = "RNA_pseudoXenium2")
vrMainFeatureType(VRBlock[["Assay2"]]) <- "RNA_pseudoXenium1"
vrMainFeatureType(VRBlock, assay = "all")

### Visium vs XeniumR1 ####

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

#### TACSTD2 ####

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
ggsave(filename = "../../../../../Documents/VoltRon Paper/Nature Methods Revision/Images/Supplementary Material/ImageAlign/images/voltron_scatterplot_TACSTD2.pdf", plot = last_plot(), 
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
# Lee's L statistic standard deviate = 110.48, p-value < 2.2e-16
# alternative hypothesis: greater
# sample estimates:
# Lee's L statistic       Expectation          Variance 
# 6.614481e-01      1.445235e-01      2.189073e-05 

#### ERBB2 ####

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
ggsave(filename = "../../../../../Documents/VoltRon Paper/Nature Methods Revision/Images/Supplementary Material/ImageAlign/images/voltron_scatterplot_ERBB2.pdf", plot = last_plot(), 
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
# Lee's L statistic standard deviate = 126.53, p-value < 2.2e-16
# alternative hypothesis: greater
# sample estimates:
# Lee's L statistic       Expectation          Variance 
# 0.7379073918      0.1461916003      0.0000218691 

### Visium vs XeniumR2 ####

# correlation plot of all genes across two assays
data1 <- vrData(VRBlock, assay = "Assay2", feat_type = "RNA")
data2 <- vrData(VRBlock, assay = "Assay2", feat_type = "RNA_pseudoXenium2")
data1 <- rowMeans(data1)
data2 <- rowMeans(data2)
datax <- data.frame(Visium = log10(data1[names(data2)]),
                    Xenium = log10(data2))
p3 <- ggplot(data = datax, mapping = aes(x = Visium, y = Xenium)) + 
  geom_point(color = "darkgreen") + 
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") + 
  geom_abline(slope = 1, intercept = 1, color = "darkgreen") +
  theme_classic() + labs(title = "Sensitivity (all panel genes)")

# correlation plot of TACSTD2 across two assays
data1 <- vrData(VRBlock, assay = "Assay2", features = "TACSTD2", feat_type = "RNA")
data2 <- vrData(VRBlock, assay = "Assay2", features = "TACSTD2", feat_type = "RNA_pseudoXenium2")
data1 <- data1[1,,drop=TRUE]
data2 <- data2[1,,drop=TRUE]
datax <- data.frame(Visium = log10(data1[names(data2)] + 1), 
                    Xenium = log10(data2 + 1))
p4 <- ggplot(data = datax, mapping = aes(x = Visium, y = Xenium)) + 
  geom_point(color = "darkgreen") + 
  theme_classic() + labs(title = "Specificity (TACSTD2)")

### visualize all ####
(p1 | p3) / (p2 | p4)
ggsave(filename = "../../allresults/BreastCancer/figures/scatterplot/voltron_XenvsVisium.pdf", 
       plot = p1 | p2 | p3 | p4,
       device = "pdf", width = 24, height = 6)

####
## compare visualization ####
####

VRBlock_assay1 <- subset(VRBlock, assay = "Assay1")
vrImages(VRBlock_assay1[["Assay1"]], name = "main_reg", channel = "H&E") <- vrImages(VRBlock, assay = "Assay2")
vrMainChannel(VRBlock_assay1[["Assay1"]]) <- "H&E"
VRBlock_assay1 <- subset(VRBlock_assay1, image = readRDS("../../allresults/BreastCancer/figures/visium_align_zoom/Visium_vis_subset.rds"))
vrSpatialPlot(VRBlock_assay1, spatial = "main_reg", channel = "H&E", color = list(`10XBlock` = "black")) + 
  theme(legend.position = "none") + 
  labs(title = "")
ggsave(filename = "../../allresults/BreastCancer/figures/visium_align_zoom/voltron_subset_breastcancervisium.pdf", 
       plot = last_plot(), width = 7, height = 7)

## visualize plot overlay ####
img1 <- vrImages(VRBlock, assay = "Assay2")
img2 <- vrImages(VRBlock, assay = "Assay1")
img1 <- image_resize(img1, geometry = magick::geometry_size_percent(20))
img2 <- image_resize(img2, geometry = magick::geometry_size_percent(20))
img2 <- image_modulate(img2, brightness = 300)
img2 <- image_negate(img2)
result <- image_composite(img1, img2, operator = "Bumpmap")
magick::image_ggplot(result)
ggsave(filename = "../../allresults/BreastCancer/figures/visium_overlay/voltron_XenvsVisium_overlay.pdf", 
       plot = last_plot(),
       device = "pdf", width = 6, height = 8)