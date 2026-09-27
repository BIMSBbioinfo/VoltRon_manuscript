library(ggplot2)
library(reshape2)
library(dplyr)

####
# Xenium vs HE ####
####

voltron_results <- read.csv("../../VoltRon/BreastCancer/results/XeniumR1_vs_HE/all_landmarks.csv")
stalign_results <- read.csv("../../STAlign/BreastCancer/results/XeniumR1_vs_HE/all_landmarks.csv")
wsireg_results <- read.csv("../../wsireg/BreastCancer/results/HE_vs_XeniumR1/all_landmarks.csv")
# all_results <- data.frame(stalign_results[,c("diff_affine", "diff_nonrigid")],
#                           diff_rigid = voltron_results[,c("diff_affine", "diff_affinetps", "diff_auto", "diff_autotps")])
# colnames(all_results) <- c("STAlign_affine", "STAlign_OT", "VoltRon_Affine", "VoltRon_Affine+TPS", "VoltRon_Auto", "VoltRon_Auto+TPS")
all_results <- data.frame(stalign = stalign_results[,c("diff_affine", "diff_nonrigid")],
                          wesireg = wsireg_results[,c("diff_affine", "diff_nl")],
                          voltron = voltron_results[,c("diff_manualaffine", "diff_affine", "diff_auto", "diff_affine_nr")])
colnames(all_results) <- c("STAlign_ManualAffine", "STAlign_OT", "WSIReg_AutoAffine", "WSIReg_NL", "VoltRon_ManualAffine", "VoltRon_AutoAffine", "VoltRon_AutoHomo", "VoltRon_AutoAffine+NR")
all_results <- reshape2::melt(all_results)
colnames(all_results) <- c("method", "diff")
all_results$diff <- log10(all_results$diff)
g1 <- ggplot(data = all_results, mapping = aes(x = method, y = diff, color = method)) + 
  geom_boxplot() + 
  geom_point(stat = "identity", position = position_jitter()) + labs(title = "Xenium (query) vs H&E (reference)") + 
  theme_classic() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

####
# Xenium vs DAPI ####
####

voltron_results <- read.csv("../../VoltRon/BreastCancer/results/XeniumR1_vs_DAPI/all_landmarks.csv")
stalign_results <- read.csv("../../STAlign/BreastCancer/results/XeniumR1_vs_DAPI/all_landmarks.csv")
wsireg_results <- read.csv("../../wsireg/BreastCancer/results/IF_vs_XeniumR1/all_landmarks.csv")
# all_results <- data.frame(stalign_results[,c("diff_affine", "diff_nonrigid")],
#                           diff_rigid = voltron_results[,c("diff_manualtps", "diff_manualhomotps", "diff_auto", "diff_autotps")])
# colnames(all_results) <- c("STAlign_affine", "STAlign_OT", "VoltRon_TPS", "VoltRon_Homo+TPS", "VoltRon_Auto", "VoltRon_Auto+TPS")
all_results <- data.frame(stalign = stalign_results[,c("diff_affine", "diff_nonrigid")],
                          wesireg = wsireg_results[,c("diff_affine", "diff_nl")],
                          voltron = voltron_results[,c("diff_manualaffine", "diff_affine", "diff_auto", "diff_affine_nr")])
colnames(all_results) <- c("STAlign_ManualAffine", "STAlign_OT", "WSIReg_AutoAffine", "WSIReg_NL", "VoltRon_ManualAffine", "VoltRon_AutoAffine", "VoltRon_AutoHomo", "VoltRon_AutoAffine+NR")
all_results <- reshape2::melt(all_results)
colnames(all_results) <- c("method", "diff")
all_results$diff <- log10(all_results$diff)
g2 <- ggplot(data = all_results, mapping = aes(x = method, y = diff, color = method)) + 
  geom_boxplot() + 
  geom_point(stat = "identity", position = position_jitter()) + labs(title = "Xenium (query) vs DAPI (reference)") + 
  theme_classic() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g2 <- g2 + ylim(0,2.7) 

# save image
ggsave(filename = "all_methods_comparison.pdf", plot = g1 | g2, width = 12, height = 7)

####
# Make images ####
####

####
## VoltRon ####
####

# import data
XeniumR1_image <- magick::image_read(paste0("../../../../../analysis/ImageAlignmentBenchmark/data/BreastCancer/Xenium_R1/Xenium/image_3.tiff"))
DAPI_image <- magick::image_read(paste0("../../../../../analysis/ImageAlignmentBenchmark/data/BreastCancer/Xenium_R1/IF/DAPI.tiff"))
HE_image <- magick::image_read(paste0("../../../../../analysis/ImageAlignmentBenchmark/data/BreastCancer/Xenium_R1/HE/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image_highres.tif"))

# voltron
ref_image <- HE_image
query_image <- XeniumR1_image
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(query_image), 
                               mapping_parameters = readRDS("../../../benchmarks/VoltRon/BreastCancer/parameters/XeniumR1_vs_HE_affine.rds"), interactive = FALSE)
img1 <- vrImages(xen_reg$registered_spat[[1]], name = "main")
img2 <- vrImages(xen_reg$registered_spat[[2]], name = "main_reg")
img1 <- image_resize(img1, geometry = magick::geometry_size_percent(20))
img2 <- image_resize(img2, geometry = magick::geometry_size_percent(20))
img2 <- image_negate(img2)
result <- image_composite(img1, img2, operator = "Bumpmap")
magick::image_ggplot(result)
ggsave(filename = "results/voltron_Xenium_vs_HE.pdf", plot = last_plot(), width = 7, height = 7)

ref_image <- DAPI_image
query_image <- XeniumR1_image
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(query_image), 
                               mapping_parameters = readRDS("../../../benchmarks/VoltRon/BreastCancer/parameters/XeniumR1_vs_DAPI_affine.rds"), interactive = FALSE)
img1 <- vrImages(xen_reg$registered_spat[[1]], name = "main")
img2 <- vrImages(xen_reg$registered_spat[[2]], name = "main_reg")
img1 <- image_resize(img1, geometry = magick::geometry_size_percent(20))
img2 <- image_resize(img2, geometry = magick::geometry_size_percent(20))
img1 <- image_modulate(img1, brightness = 300)
img2 <- image_modulate(img2, brightness = 200)
result <- image_composite(img1, img2, operator = "Blend")
magick::image_ggplot(result)
ggsave(filename = "results/voltron_Xenium_vs_DAPI.pdf", plot = last_plot(), width = 7, height = 7)

####
## WSIReg ####
####

img1 <- HE_image
img2 <- RBioFormats::read.image("../../../benchmarks/wsireg/BreastCancer/results/HE_vs_XeniumR1/my_reg_project-modality_fluo_to_modality_brightfield_registered.ome.tiff", 
                               series = 1, resolution = 1)
img2 <- magick::image_read(as.raster(img2))
img1 <- image_resize(img1, geometry = magick::geometry_size_percent(20))
img2 <- image_resize(img2, geometry = magick::geometry_size_percent(20))
img2 <- image_negate(img2)
result <- image_composite(img1, img2, operator = "Bumpmap")
magick::image_ggplot(result)
ggsave(filename = "results/wsireg_Xenium_vs_HE.pdf", plot = last_plot(), width = 7, height = 7)

img1 <- RBioFormats::read.image("../../../benchmarks/wsireg/BreastCancer/results/IF_vs_XeniumR1/my_reg_project-modality_brightfield_registered.ome.tiff", 
                                series = 1, resolution = 1)
img1 <- magick::image_read(as.raster(img1))
img2 <- RBioFormats::read.image("../../../benchmarks/wsireg/BreastCancer/results/IF_vs_XeniumR1/my_reg_project-modality_fluo_to_modality_brightfield_registered.ome.tiff", 
                                series = 1, resolution = 1)
img2 <- magick::image_read(as.raster(img2))
img1 <- image_resize(img1, geometry = magick::geometry_size_percent(20))
img2 <- image_resize(img2, geometry = magick::geometry_size_percent(20))
img1 <- image_modulate(img1, brightness = 300)
img2 <- image_modulate(img2, brightness = 200)
result <- image_composite(img1, img2, operator = "Blend")
magick::image_ggplot(result)
ggsave(filename = "results/wsireg_Xenium_vs_DAPI.pdf", plot = last_plot(), width = 7, height = 7)

####
## STAlign ####
####

img1 <- HE_image
points <- read.csv("../../../benchmarks/STAlign/BreastCancer/results/XeniumR1_vs_HE/transformed_points.csv")
points <- points[,c(3,2)]
colnames(points) <- c("x", "y")
points$y <- magick::image_info(img1)$height - points$y
img1 <- image_resize(img1, geometry = magick::geometry_size_percent(20))
points <- points*0.20
magick::image_ggplot(img1) + 
  geom_point(data = points, mapping = aes(x = x, y = y), color = "black", size = 0.05, alpha = 0.5)
ggsave(filename = "results/stalign_Xenium_vs_HE.png", device = "png", plot = last_plot(), width = 7, height = 7)

img1 <- DAPI_image
points <- read.csv("../../../benchmarks/STAlign/BreastCancer/results/XeniumR1_vs_DAPI/transformed_points.csv")
points <- points[,c(3,2)]
colnames(points) <- c("x", "y")
points$y <- magick::image_info(img1)$height - points$y
img1 <- image_resize(img1, geometry = magick::geometry_size_percent(20))
img1 <- image_modulate(img1, brightness = 300)
points <- points*0.20
g1 <- magick::image_ggplot(img1) + 
  stat_bin_2d(data = points, mapping = aes(x = x, y = y), color = "red", size = 0.05, alpha = 0.5, bins = 200)
datax <- g1@layers$stat_bin_2d$data
magick::image_ggplot(img1) + 
  geom_point(data = datax, mapping = aes(x = x, y = y), color = "red", size = 0.05, alpha = 0.5)
ggsave(filename = "results/stalign_Xenium_vs_DAPI.png", device = "png", plot = last_plot(), width = 7, height = 7)