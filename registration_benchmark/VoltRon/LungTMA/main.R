# library(VoltRon)
library(ggpubr)
library(magick)

####
# Get TMA ####
####

nimages <- 8
DAPI_image_list <- list(
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/Xenium/lungTMA_Assay1.tif"),
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/Xenium/lungTMA_Assay2.tif"),
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/Xenium/lungTMA_Assay3.tif"),
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/Xenium/lungTMA_Assay4.tif"),
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/Xenium/lungTMA_Assay5.tif"),
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/Xenium/lungTMA_Assay6.tif"),
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/Xenium/lungTMA_Assay7.tif"),
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/Xenium/lungTMA_Assay8.tif")
)

HE_image_list <- list(
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/HE/lungTMA_1.jpg"),
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/HE/lungTMA_2.jpg"),
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/HE/lungTMA_3.jpg"),
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/HE/lungTMA_4.jpg"),
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/HE/lungTMA_5.jpg"),
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/HE/lungTMA_6.jpg"),
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/HE/lungTMA_7.jpg"),
  magick::image_read("../../../../ImageAlignmentBenchmark/data/LungTMA/HE/lungTMA_8.jpg")
)

# get image resolutions
cat(matrix(sapply(DAPI_image_list, function(x) paste(image_info(x)[,c("width", "height")][1,], collapse = "x"), simplify = TRUE), ncol = 1))
cat(matrix(sapply(HE_image_list, function(x) paste(image_info(x)[,c("width", "height")][1,], collapse = "x"), simplify = TRUE), ncol = 1))

####
# Lung TMA 1 ####
####

####
## get parameters ####
####

ref_image <- HE_image_list[[1]]
query_image <- DAPI_image_list[[1]]
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 400)), 
                               mapping_parameters = readRDS("parameters/lungTMA1_HE_vs_DAPI_affine_2.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA1_HE_vs_DAPI_affine_2.rds")

# xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
#                                query_spatdata = importImageData(image_modulate(query_image, brightness = 400)), 
#                                mapping_parameters = readRDS("parameters/lungTMA1_HE_vs_DAPI_homography.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA1_HE_vs_DAPI_homography.rds")

####
## transformation of validation landmarks ####
####

# landmarks and images
ref_image <- HE_image_list[[1]]
query_image <- magick::image_rotate(image_modulate(DAPI_image_list[[1]], brightness = 400), degrees = 90)
HE_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA1/HE_landmarks_corrected.csv", row.names = 1) 
Xenium_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA1/Xenium_landmarks_corrected.csv", row.names = 1) 
imageinfo_query <- magick::image_info(query_image)
imageinfo_ref <- magick::image_info(ref_image)

# transform validation markers (affine)
mapping_parameters <- readRDS("parameters/lungTMA1_HE_vs_DAPI_affine_2.rds")
Xenium_landmarks_affine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine$y <- imageinfo_query$height - Xenium_landmarks_affine$y
Xenium_landmarks_affine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine$y <- imageinfo_ref$height - Xenium_landmarks_affine$y
diff_affine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine[,c("x", "y")])^2))

# transform validation markers (affine + BSpline)
mapping <- readRDS("parameters/lungTMA1_HE_vs_DAPI_affine_2.rds")
mapping$Method <- "Affine + Non-Rigid"
mapping$nonrigid <- "BSpline (SimpleITK)"
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(DAPI_image_list[[1]], brightness = 400)),
                               mapping_parameters = mapping, interactive = TRUE)
mapping_parameters <- xen_reg$mapping_parameters
Xenium_landmarks_affine_nr <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine_nr$y <- imageinfo_query$height - Xenium_landmarks_affine_nr$y
Xenium_landmarks_affine_nr[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine_nr), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine_nr$y <- imageinfo_ref$height - Xenium_landmarks_affine_nr$y
diff_affine_nr <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine_nr[,c("x", "y")])^2))

# compare markers
results <- data.frame(HE_landmarks[,c("x","y")], diff_affine, diff_affine, diff_affine, diff_affine_nr)
colnames(results) <- c('x_landmark', 'y_landmark', "diff_manualaffine", "diff_affine", "diff_auto", "diff_affine_nr")
dir.create("results/lungTMA1/")
write.csv(results, file = "results/lungTMA1/all_landmarks.csv", row.names = FALSE, quote = FALSE)

####
# Lung TMA 2 ####
####

####
## get parameters ####
####

ref_image <- HE_image_list[[2]]
query_image <- DAPI_image_list[[2]]
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)),
                               mapping_parameters = readRDS("parameters/lungTMA2_HE_vs_DAPI_homography.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA2_HE_vs_DAPI_homography.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)),
                               mapping_parameters = readRDS("parameters/lungTMA2_HE_vs_DAPI_affine.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA2_HE_vs_DAPI_affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)),
                               mapping_parameters = readRDS("parameters/lungTMA2_HE_vs_DAPI_affinetps.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA2_HE_vs_DAPI_affinetps.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA2_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA2_HE_vs_DAPI_manual_Affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA2_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA2_HE_vs_DAPI_manual_TPS.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA2_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA2_HE_vs_DAPI_manual_Affine_TPS.rds")

####
## transformation of validation landmarks ####
####

# landmarks and images
ref_image <- HE_image_list[[2]]
query_image <- magick::image_rotate(image_modulate(DAPI_image_list[[2]], brightness = 200), degrees = 90)
HE_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA2/HE_landmarks_corrected.csv", row.names = 1) 
Xenium_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA2/Xenium_landmarks_corrected.csv", row.names = 1) 
imageinfo_query <- magick::image_info(query_image)
imageinfo_ref <- magick::image_info(ref_image)

# transform validation markers (auto)
mapping_parameters <- readRDS("parameters/lungTMA2_HE_vs_DAPI_homography.rds")
Xenium_landmarks_auto <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_auto$y <- imageinfo_query$height - Xenium_landmarks_auto$y
Xenium_landmarks_auto[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_auto), mapping_parameters$mapping$`2`)
Xenium_landmarks_auto$y <- imageinfo_ref$height - Xenium_landmarks_auto$y
diff_auto <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_auto[,c("x", "y")])^2))

# transform validation markers (affine)
mapping_parameters <- readRDS("parameters/lungTMA2_HE_vs_DAPI_affine.rds")
Xenium_landmarks_affine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine$y <- imageinfo_query$height - Xenium_landmarks_affine$y
Xenium_landmarks_affine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine$y <- imageinfo_ref$height - Xenium_landmarks_affine$y
diff_affine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine[,c("x", "y")])^2))

# transform validation markers (affine + TPS)
mapping_parameters <- readRDS("parameters/lungTMA2_HE_vs_DAPI_affinetps.rds")
Xenium_landmarks_affinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affinetps$y <- imageinfo_query$height - Xenium_landmarks_affinetps$y
Xenium_landmarks_affinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_affinetps$y <- imageinfo_ref$height - Xenium_landmarks_affinetps$y
diff_affinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affinetps[,c("x", "y")])^2))

# transform validation markers (manual homo + tps)
mapping_parameters <- readRDS("parameters/lungTMA2_HE_vs_DAPI_manual_Affine.rds")
Xenium_landmarks_manualaffine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffine$y <- imageinfo_query$height - Xenium_landmarks_manualaffine$y
Xenium_landmarks_manualaffine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffine), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffine$y <- imageinfo_ref$height - Xenium_landmarks_manualaffine$y
diff_manualaffine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffine[,c("x", "y")])^2))

# transform validation markers (manual tps)
mapping_parameters <- readRDS("parameters/lungTMA2_HE_vs_DAPI_manual_TPS.rds")
Xenium_landmarks_manualtps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualtps$y <- imageinfo_query$height - Xenium_landmarks_manualtps$y
Xenium_landmarks_manualtps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualtps), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualtps$y <- imageinfo_ref$height - Xenium_landmarks_manualtps$y
diff_manualtps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualtps[,c("x", "y")])^2))

# transform validation markers (manual homo + tps)
mapping_parameters <- readRDS("parameters/lungTMA2_HE_vs_DAPI_manual_Affine_TPS.rds")
Xenium_landmarks_manualaffinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffinetps$y <- imageinfo_query$height - Xenium_landmarks_manualaffinetps$y
Xenium_landmarks_manualaffinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffinetps$y <- imageinfo_ref$height - Xenium_landmarks_manualaffinetps$y
diff_manualaffinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffinetps[,c("x", "y")])^2))

# transform validation markers (affine + BSpline)
mapping <- readRDS("parameters/lungTMA2_HE_vs_DAPI_affine.rds")
mapping$Method <- "Affine + Non-Rigid"
mapping$nonrigid <- "BSpline (SimpleITK)"
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(DAPI_image_list[[2]], brightness = 200)),
                               mapping_parameters = mapping, interactive = TRUE)
mapping_parameters <- xen_reg$mapping_parameters
Xenium_landmarks_affine_nr <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine_nr$y <- imageinfo_query$height - Xenium_landmarks_affine_nr$y
Xenium_landmarks_affine_nr[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine_nr), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine_nr$y <- imageinfo_ref$height - Xenium_landmarks_affine_nr$y
diff_affine_nr <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine_nr[,c("x", "y")])^2))

# compare markers
results <- data.frame(HE_landmarks[,c("x","y")], diff_manualaffine, diff_manualtps, diff_manualaffinetps, diff_affine, diff_affinetps, diff_auto, diff_affine_nr)
colnames(results) <- c('x_landmark', 'y_landmark', "diff_manualaffine", "diff_manualtps", "diff_manualaffinetps", "diff_affine", "diff_affinetps", "diff_auto", "diff_affine_nr")
dir.create("results/lungTMA2/")
write.csv(results, file = "results/lungTMA2/all_landmarks.csv", row.names = FALSE, quote = FALSE)

####
# Lung TMA 3 ####
####

####
## get parameters ####
####

ref_image <- HE_image_list[[3]]
query_image <- DAPI_image_list[[3]]
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA3_HE_vs_DAPI_homography.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA3_HE_vs_DAPI_homography.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA3_HE_vs_DAPI_affine.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA3_HE_vs_DAPI_affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA3_HE_vs_DAPI_affinetps.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA3_HE_vs_DAPI_affinetps.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA3_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA3_HE_vs_DAPI_manual_Affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA3_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA3_HE_vs_DAPI_manual_TPS.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA3_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA3_HE_vs_DAPI_manual_Affine_TPS.rds")

####
## image transformation ####
####

# apply transformation to image homography
mapping_parameters <- readRDS("parameters/lungTMA3_HE_vs_DAPI_homography.rds")
image3 <- getRcppWarpImage(ref_image = ref_image, 
                           query_image = query_image,
                           mapping = mapping_parameters$mapping$`2`)
image_view_list <- list(rep(magick::image_resize(HE_image_list[[3]], geometry = "5000x"),1),
                        rep(magick::image_resize(image3, geometry = "5000x"),1))
image_view_list %>% magick::image_join()

####
## transformation of validation landmarks ####
####

# landmarks and images
ref_image <- HE_image_list[[3]]
query_image <- magick::image_rotate(image_modulate(DAPI_image_list[[3]], brightness = 200), degrees = 90)
HE_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA3/HE_landmarks_corrected.csv", row.names = 1) 
Xenium_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA3/Xenium_landmarks_corrected.csv", row.names = 1) 
imageinfo_query <- magick::image_info(query_image)
imageinfo_ref <- magick::image_info(ref_image)

# transform validation markers (auto)
mapping_parameters <- readRDS("parameters/lungTMA3_HE_vs_DAPI_homography.rds")
Xenium_landmarks_auto <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_auto$y <- imageinfo_query$height - Xenium_landmarks_auto$y
Xenium_landmarks_auto[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_auto), mapping_parameters$mapping$`2`)
Xenium_landmarks_auto$y <- imageinfo_ref$height - Xenium_landmarks_auto$y
diff_auto <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_auto[,c("x", "y")])^2))

# transform validation markers (affine)
mapping_parameters <- readRDS("parameters/lungTMA3_HE_vs_DAPI_affine.rds")
Xenium_landmarks_affine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine$y <- imageinfo_query$height - Xenium_landmarks_affine$y
Xenium_landmarks_affine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine$y <- imageinfo_ref$height - Xenium_landmarks_affine$y
diff_affine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine[,c("x", "y")])^2))

# transform validation markers (affine + TPS)
mapping_parameters <- readRDS("parameters/lungTMA3_HE_vs_DAPI_affinetps.rds")
Xenium_landmarks_affinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affinetps$y <- imageinfo_query$height - Xenium_landmarks_affinetps$y
Xenium_landmarks_affinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_affinetps$y <- imageinfo_ref$height - Xenium_landmarks_affinetps$y
diff_affinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affinetps[,c("x", "y")])^2))

# transform validation markers (manual homo + tps)
mapping_parameters <- readRDS("parameters/lungTMA3_HE_vs_DAPI_manual_Affine.rds")
Xenium_landmarks_manualaffine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffine$y <- imageinfo_query$height - Xenium_landmarks_manualaffine$y
Xenium_landmarks_manualaffine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffine), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffine$y <- imageinfo_ref$height - Xenium_landmarks_manualaffine$y
diff_manualaffine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffine[,c("x", "y")])^2))

# transform validation markers (manual tps)
mapping_parameters <- readRDS("parameters/lungTMA3_HE_vs_DAPI_manual_TPS.rds")
Xenium_landmarks_manualtps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualtps$y <- imageinfo_query$height - Xenium_landmarks_manualtps$y
Xenium_landmarks_manualtps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualtps), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualtps$y <- imageinfo_ref$height - Xenium_landmarks_manualtps$y
diff_manualtps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualtps[,c("x", "y")])^2))

# transform validation markers (manual homo + tps)
mapping_parameters <- readRDS("parameters/lungTMA3_HE_vs_DAPI_manual_Affine_TPS.rds")
Xenium_landmarks_manualaffinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffinetps$y <- imageinfo_query$height - Xenium_landmarks_manualaffinetps$y
Xenium_landmarks_manualaffinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffinetps$y <- imageinfo_ref$height - Xenium_landmarks_manualaffinetps$y
diff_manualaffinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffinetps[,c("x", "y")])^2))

# transform validation markers (affine + BSpline)
mapping <- readRDS("parameters/lungTMA3_HE_vs_DAPI_affine.rds")
mapping$Method <- "Affine + Non-Rigid"
mapping$nonrigid <- "BSpline (SimpleITK)"
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(DAPI_image_list[[3]], brightness = 200)),
                               mapping_parameters = mapping, interactive = TRUE)
mapping_parameters <- xen_reg$mapping_parameters
Xenium_landmarks_affine_nr <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine_nr$y <- imageinfo_query$height - Xenium_landmarks_affine_nr$y
Xenium_landmarks_affine_nr[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine_nr), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine_nr$y <- imageinfo_ref$height - Xenium_landmarks_affine_nr$y
diff_affine_nr <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine_nr[,c("x", "y")])^2))

# compare markers
results <- data.frame(HE_landmarks[,c("x","y")], diff_manualaffine, diff_manualtps, diff_manualaffinetps, diff_affine, diff_affinetps, diff_auto, diff_affine_nr)
colnames(results) <- c('x_landmark', 'y_landmark', "diff_manualaffine", "diff_manualtps", "diff_manualaffinetps", "diff_affine", "diff_affinetps", "diff_auto", "diff_affine_nr")
dir.create("results/lungTMA3/")
write.csv(results, file = "results/lungTMA3/all_landmarks.csv", row.names = FALSE, quote = FALSE)
  
####
# Lung TMA 4 ####
####

####
## get parameters ####
####

ref_image <- HE_image_list[[4]]
query_image <- DAPI_image_list[[4]]
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA4_HE_vs_DAPI_homography.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA4_HE_vs_DAPI_homography.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA4_HE_vs_DAPI_affine.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA4_HE_vs_DAPI_affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA4_HE_vs_DAPI_affinetps.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA4_HE_vs_DAPI_affinetps.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA4_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA4_HE_vs_DAPI_manual_Affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA4_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA4_HE_vs_DAPI_manual_TPS.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA4_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA4_HE_vs_DAPI_manual_Affine_TPS.rds")

####
## image transformation ####
####

tmp <- readRDS("parameters/lungTMA4_HE_vs_DAPI_affinetps.rds")
keypoints <- tmp$keypoints
refpoints <- tmp$mapping$`2`[[1]][[2]][[1]]
refpoints <- as_tibble(data.frame(KeyPoint = 1:nrow(refpoints), x = refpoints[,1], y = refpoints[,2]))
keypoints$`1-2`$ref <- refpoints
querypoints <- tmp$mapping$`2`[[1]][[2]][[2]]
querypoints <- as_tibble(data.frame(KeyPoint = 1:nrow(querypoints), x = querypoints[,1], y = querypoints[,2]))
keypoints$`1-2`$query <- querypoints
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(query_image, degrees = 90)), 
                               keypoints = keypoints)

# apply transformation to image homography
mapping_parameters <- readRDS("parameters/lungTMA4_HE_vs_DAPI_affinetps.rds")
image4 <- getRcppWarpImage(ref_image = ref_image, 
                           query_image = image_rotate(query_image, degrees = 90),
                           mapping = mapping_parameters$mapping$`2`)
image_view_list <- list(rep(magick::image_resize(HE_image_list[[4]], geometry = "2000x"),1),
                        rep(magick::image_resize(image4, geometry = "2000x"),1))
image_view_list %>% magick::image_join()

# apply transformation to sift keypoints
tmp <- readRDS("parameters/lungTMA4_HE_vs_DAPI_affinetps.rds")
refpoints <- tmp$mapping$`2`[[1]][[2]][[1]]
refpoints <- data.frame(refpoints)
colnames(refpoints) <- c("x","y")
querypoints <- tmp$mapping$`2`[[1]][[2]][[2]]
querypoints_reg <- applyMapping(querypoints, mapping = tmp$mapping$`2`)
querypoints_reg <- data.frame(querypoints_reg)
colnames(querypoints_reg) <- c("x","y")
g1 <- magick::image_ggplot(magick::image_resize(HE_image_list[[4]], geometry = magick::geometry_size_percent(10))) + 
  geom_point(mapping = aes(x = x, y = y), data = refpoints*0.1, color = "yellow") + 
  geom_point(mapping = aes(x = x, y = y), data = querypoints_reg*0.1, color = "purple") 

####
## transformation of validation landmarks ####
####

# landmarks and images
ref_image <- HE_image_list[[4]]
query_image <- magick::image_rotate(image_modulate(DAPI_image_list[[4]], brightness = 200), degrees = 90)
HE_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA4/HE_landmarks_corrected.csv", row.names = 1) 
Xenium_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA4/Xenium_landmarks_corrected.csv", row.names = 1) 
imageinfo_query <- magick::image_info(query_image)
imageinfo_ref <- magick::image_info(ref_image)

# transform validation markers (auto)
mapping_parameters <- readRDS("parameters/lungTMA4_HE_vs_DAPI_homography.rds")
Xenium_landmarks_auto <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_auto$y <- imageinfo_query$height - Xenium_landmarks_auto$y
Xenium_landmarks_auto[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_auto), mapping_parameters$mapping$`2`)
Xenium_landmarks_auto$y <- imageinfo_ref$height - Xenium_landmarks_auto$y
diff_auto <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_auto[,c("x", "y")])^2))

# transform validation markers (affine)
mapping_parameters <- readRDS("parameters/lungTMA4_HE_vs_DAPI_affine.rds")
Xenium_landmarks_affine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine$y <- imageinfo_query$height - Xenium_landmarks_affine$y
Xenium_landmarks_affine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine$y <- imageinfo_ref$height - Xenium_landmarks_affine$y
diff_affine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine[,c("x", "y")])^2))

# transform validation markers (affine + TPS)
mapping_parameters <- readRDS("parameters/lungTMA4_HE_vs_DAPI_affinetps.rds")
Xenium_landmarks_affinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affinetps$y <- imageinfo_query$height - Xenium_landmarks_affinetps$y
Xenium_landmarks_affinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_affinetps$y <- imageinfo_ref$height - Xenium_landmarks_affinetps$y
diff_affinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affinetps[,c("x", "y")])^2))

# transform validation markers (manual affine)
mapping_parameters <- readRDS("parameters/lungTMA4_HE_vs_DAPI_manual_Affine.rds")
Xenium_landmarks_manualaffine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffine$y <- imageinfo_query$height - Xenium_landmarks_manualaffine$y
Xenium_landmarks_manualaffine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffine), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffine$y <- imageinfo_ref$height - Xenium_landmarks_manualaffine$y
diff_manualaffine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffine[,c("x", "y")])^2))

# transform validation markers (manual tps)
mapping_parameters <- readRDS("parameters/lungTMA4_HE_vs_DAPI_manual_TPS.rds")
Xenium_landmarks_manualtps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualtps$y <- imageinfo_query$height - Xenium_landmarks_manualtps$y
Xenium_landmarks_manualtps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualtps), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualtps$y <- imageinfo_ref$height - Xenium_landmarks_manualtps$y
diff_manualtps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualtps[,c("x", "y")])^2))

# transform validation markers (manual homo + tps)
mapping_parameters <- readRDS("parameters/lungTMA4_HE_vs_DAPI_manual_Affine_TPS.rds")
Xenium_landmarks_manualaffinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffinetps$y <- imageinfo_query$height - Xenium_landmarks_manualaffinetps$y
Xenium_landmarks_manualaffinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffinetps$y <- imageinfo_ref$height - Xenium_landmarks_manualaffinetps$y
diff_manualaffinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffinetps[,c("x", "y")])^2))

# transform validation markers (affine + BSpline)
mapping_parameters <- readRDS("parameters/lungTMA4_HE_vs_DAPI_affine.rds")
mapping$Method <- "Affine + Non-Rigid"
mapping$nonrigid <- "BSpline (SimpleITK)"
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(DAPI_image_list[[4]], brightness = 200)),
                               mapping_parameters = mapping, interactive = TRUE)
mapping_parameters <- xen_reg$mapping_parameters
Xenium_landmarks_affine_nr <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine_nr$y <- imageinfo_query$height - Xenium_landmarks_affine_nr$y
Xenium_landmarks_affine_nr[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine_nr), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine_nr$y <- imageinfo_ref$height - Xenium_landmarks_affine_nr$y
diff_affine_nr <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine_nr[,c("x", "y")])^2))

# compare markers
results <- data.frame(HE_landmarks[,c("x","y")], diff_manualaffine, diff_manualtps, diff_manualaffinetps, diff_affine, diff_affinetps, diff_auto, diff_affine_nr)
colnames(results) <- c('x_landmark', 'y_landmark', "diff_manualaffine", "diff_manualtps", "diff_manualaffinetps", "diff_affine", "diff_affinetps", "diff_auto", "diff_affine_nr")
dir.create("results/lungTMA4/")
write.csv(results, file = "results/lungTMA4/all_landmarks.csv", row.names = FALSE, quote = FALSE)
  
####
# Lung TMA 5 ####
####

####
## get parameters ####
####

ref_image <- HE_image_list[[5]]
query_image <- DAPI_image_list[[5]]
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS(file = "parameters/lungTMA5_HE_vs_DAPI_homography.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA5_HE_vs_DAPI_homography.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS(file = "parameters/lungTMA5_HE_vs_DAPI_affine.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA5_HE_vs_DAPI_affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS(file = "parameters/lungTMA5_HE_vs_DAPI_affinetps.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA5_HE_vs_DAPI_affinetps.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA5_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA5_HE_vs_DAPI_manual_Affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA5_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA5_HE_vs_DAPI_manual_TPS.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA5_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA5_HE_vs_DAPI_manual_Affine_TPS.rds")

####
## transformation of validation landmarks ####
####

# landmarks and images
ref_image <- HE_image_list[[5]]
query_image <- magick::image_rotate(image_modulate(DAPI_image_list[[5]], brightness = 200), degrees = 90)
HE_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA5/HE_landmarks_corrected.csv", row.names = 1) 
Xenium_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA5/Xenium_landmarks_corrected.csv", row.names = 1) 
imageinfo_query <- magick::image_info(query_image)
imageinfo_ref <- magick::image_info(ref_image)

# transform validation markers (auto)
mapping_parameters <- readRDS("parameters/lungTMA5_HE_vs_DAPI_homography.rds")
Xenium_landmarks_auto <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_auto$y <- imageinfo_query$height - Xenium_landmarks_auto$y
Xenium_landmarks_auto[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_auto), mapping_parameters$mapping$`2`)
Xenium_landmarks_auto$y <- imageinfo_ref$height - Xenium_landmarks_auto$y
diff_auto <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_auto[,c("x", "y")])^2))

# transform validation markers (affine)
mapping_parameters <- readRDS("parameters/lungTMA5_HE_vs_DAPI_affine.rds")
Xenium_landmarks_affine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine$y <- imageinfo_query$height - Xenium_landmarks_affine$y
Xenium_landmarks_affine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine$y <- imageinfo_ref$height - Xenium_landmarks_affine$y
diff_affine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine[,c("x", "y")])^2))

# transform validation markers (affine + TPS)
mapping_parameters <- readRDS("parameters/lungTMA5_HE_vs_DAPI_affinetps.rds")
Xenium_landmarks_affinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affinetps$y <- imageinfo_query$height - Xenium_landmarks_affinetps$y
Xenium_landmarks_affinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_affinetps$y <- imageinfo_ref$height - Xenium_landmarks_affinetps$y
diff_affinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affinetps[,c("x", "y")])^2))

# transform validation markers (manual homo + tps)
mapping_parameters <- readRDS("parameters/lungTMA5_HE_vs_DAPI_manual_Affine.rds")
Xenium_landmarks_manualaffine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffine$y <- imageinfo_query$height - Xenium_landmarks_manualaffine$y
Xenium_landmarks_manualaffine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffine), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffine$y <- imageinfo_ref$height - Xenium_landmarks_manualaffine$y
diff_manualaffine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffine[,c("x", "y")])^2))

# transform validation markers (manual tps)
mapping_parameters <- readRDS("parameters/lungTMA5_HE_vs_DAPI_manual_TPS.rds")
Xenium_landmarks_manualtps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualtps$y <- imageinfo_query$height - Xenium_landmarks_manualtps$y
Xenium_landmarks_manualtps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualtps), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualtps$y <- imageinfo_ref$height - Xenium_landmarks_manualtps$y
diff_manualtps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualtps[,c("x", "y")])^2))

# transform validation markers (manual homo + tps)
mapping_parameters <- readRDS("parameters/lungTMA5_HE_vs_DAPI_manual_Affine_TPS.rds")
Xenium_landmarks_manualaffinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffinetps$y <- imageinfo_query$height - Xenium_landmarks_manualaffinetps$y
Xenium_landmarks_manualaffinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffinetps$y <- imageinfo_ref$height - Xenium_landmarks_manualaffinetps$y
diff_manualaffinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffinetps[,c("x", "y")])^2))

# transform validation markers (affine + BSpline)
mapping <- readRDS("parameters/lungTMA5_HE_vs_DAPI_affine.rds")
mapping$Method <- "Affine + Non-Rigid"
mapping$nonrigid <- "BSpline (SimpleITK)"
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(DAPI_image_list[[5]], brightness = 200)),
                               mapping_parameters = mapping, interactive = TRUE)
mapping_parameters <- xen_reg$mapping_parameters
Xenium_landmarks_affine_nr <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine_nr$y <- imageinfo_query$height - Xenium_landmarks_affine_nr$y
Xenium_landmarks_affine_nr[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine_nr), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine_nr$y <- imageinfo_ref$height - Xenium_landmarks_affine_nr$y
diff_affine_nr <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine_nr[,c("x", "y")])^2))

# compare markers
results <- data.frame(HE_landmarks[,c("x","y")], diff_manualaffine, diff_manualtps, diff_manualaffinetps, diff_affine, diff_affinetps, diff_auto, diff_affine_nr)
colnames(results) <- c('x_landmark', 'y_landmark', "diff_manualaffine", "diff_manualtps", "diff_manualaffinetps", "diff_affine", "diff_affinetps", "diff_auto", "diff_affine_nr")
dir.create("results/lungTMA5/")
write.csv(results, file = "results/lungTMA5/all_landmarks.csv", row.names = FALSE, quote = FALSE)

####
# Lung TMA 6 ####
####

####
## get parameters ####
####

ref_image <- HE_image_list[[6]]
query_image <- DAPI_image_list[[6]]
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA6_HE_vs_DAPI_homography.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA6_HE_vs_DAPI_homography.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA6_HE_vs_DAPI_affine.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA6_HE_vs_DAPI_affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA6_HE_vs_DAPI_affinetps.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA6_HE_vs_DAPI_affinetps.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA6_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA6_HE_vs_DAPI_manual_Affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA6_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA6_HE_vs_DAPI_manual_TPS.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA6_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA6_HE_vs_DAPI_manual_Affine_TPS.rds")

####
## image transformation ####
####

# apply transformation to image homography
mapping_parameters <- readRDS("parameters/lungTMA6_HE_vs_DAPI_homography.rds")
image4 <- getRcppWarpImage(ref_image = ref_image, 
                           query_image = query_image,
                           mapping = mapping_parameters$mapping$`2`)
image_view_list <- list(rep(magick::image_resize(HE_image_list[[4]], geometry = "5000x"),1),
                        rep(magick::image_resize(image4, geometry = "5000x"),1))
image_view_list %>% magick::image_join()

####
## transformation of validation landmarks ####
####

# landmarks and images
ref_image <- HE_image_list[[6]]
query_image <- magick::image_rotate(image_modulate(DAPI_image_list[[6]], brightness = 200), degrees = 90)
HE_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA6/HE_landmarks_corrected.csv", row.names = 1) 
Xenium_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA6/Xenium_landmarks_corrected.csv", row.names = 1) 
imageinfo_query <- magick::image_info(query_image)
imageinfo_ref <- magick::image_info(ref_image)

# transform validation markers (auto)
mapping_parameters <- readRDS("parameters/lungTMA6_HE_vs_DAPI_homography.rds")
Xenium_landmarks_auto <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_auto$y <- imageinfo_query$height - Xenium_landmarks_auto$y
Xenium_landmarks_auto[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_auto), mapping_parameters$mapping$`2`)
Xenium_landmarks_auto$y <- imageinfo_ref$height - Xenium_landmarks_auto$y
diff_auto <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_auto[,c("x", "y")])^2))

# transform validation markers (affine)
mapping_parameters <- readRDS("parameters/lungTMA6_HE_vs_DAPI_affine.rds")
Xenium_landmarks_affine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine$y <- imageinfo_query$height - Xenium_landmarks_affine$y
Xenium_landmarks_affine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine$y <- imageinfo_ref$height - Xenium_landmarks_affine$y
diff_affine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine[,c("x", "y")])^2))

# transform validation markers (affine + TPS)
mapping_parameters <- readRDS("parameters/lungTMA6_HE_vs_DAPI_affinetps.rds")
Xenium_landmarks_affinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affinetps$y <- imageinfo_query$height - Xenium_landmarks_affinetps$y
Xenium_landmarks_affinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_affinetps$y <- imageinfo_ref$height - Xenium_landmarks_affinetps$y
diff_affinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affinetps[,c("x", "y")])^2))

# transform validation markers (manual homo + tps)
mapping_parameters <- readRDS("parameters/lungTMA6_HE_vs_DAPI_manual_Affine.rds")
Xenium_landmarks_manualaffine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffine$y <- imageinfo_query$height - Xenium_landmarks_manualaffine$y
Xenium_landmarks_manualaffine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffine), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffine$y <- imageinfo_ref$height - Xenium_landmarks_manualaffine$y
diff_manualaffine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffine[,c("x", "y")])^2))

# transform validation markers (manual tps)
mapping_parameters <- readRDS("parameters/lungTMA6_HE_vs_DAPI_manual_TPS.rds")
Xenium_landmarks_manualtps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualtps$y <- imageinfo_query$height - Xenium_landmarks_manualtps$y
Xenium_landmarks_manualtps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualtps), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualtps$y <- imageinfo_ref$height - Xenium_landmarks_manualtps$y
diff_manualtps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualtps[,c("x", "y")])^2))

# transform validation markers (manual affine + tps)
mapping_parameters <- readRDS("parameters/lungTMA6_HE_vs_DAPI_manual_Affine_TPS.rds")
Xenium_landmarks_manualaffinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffinetps$y <- imageinfo_query$height - Xenium_landmarks_manualaffinetps$y
Xenium_landmarks_manualaffinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffinetps$y <- imageinfo_ref$height - Xenium_landmarks_manualaffinetps$y
diff_manualaffinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffinetps[,c("x", "y")])^2))

# transform validation markers (affine + BSpline)
mapping <- readRDS("parameters/lungTMA6_HE_vs_DAPI_affine.rds")
mapping$Method <- "Affine + Non-Rigid"
mapping$nonrigid <- "BSpline (SimpleITK)"
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(DAPI_image_list[[6]], brightness = 200)),
                               mapping_parameters = mapping, interactive = TRUE)
mapping_parameters <- xen_reg$mapping_parameters
Xenium_landmarks_affine_nr <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine_nr$y <- imageinfo_query$height - Xenium_landmarks_affine_nr$y
Xenium_landmarks_affine_nr[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine_nr), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine_nr$y <- imageinfo_ref$height - Xenium_landmarks_affine_nr$y
diff_affine_nr <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine_nr[,c("x", "y")])^2))

# compare markers
results <- data.frame(HE_landmarks[,c("x","y")], diff_manualaffine, diff_manualtps, diff_manualaffinetps, diff_affine, diff_affinetps, diff_auto, diff_affine_nr)
colnames(results) <- c('x_landmark', 'y_landmark', "diff_manualaffine", "diff_manualtps", "diff_manualaffinetps", "diff_affine", "diff_affinetps", "diff_auto", "diff_affine_nr")
dir.create("results/lungTMA6/")
write.csv(results, file = "results/lungTMA6/all_landmarks.csv", row.names = FALSE, quote = FALSE)

####
# Lung TMA 7 ####
####

####
## get parameters ####
####

ref_image <- HE_image_list[[7]]
query_image <- DAPI_image_list[[7]]
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA7_HE_vs_DAPI_homography.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA7_HE_vs_DAPI_homography.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA7_HE_vs_DAPI_affine.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA7_HE_vs_DAPI_affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA7_HE_vs_DAPI_affinetps.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA7_HE_vs_DAPI_affinetps.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA7_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA7_HE_vs_DAPI_manual_Affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA7_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA7_HE_vs_DAPI_manual_TPS.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA7_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA7_HE_vs_DAPI_manual_Affine_TPS.rds")

####
## transformation of validation landmarks ####
####

# landmarks and images
ref_image <- HE_image_list[[7]]
query_image <- magick::image_rotate(image_modulate(DAPI_image_list[[7]], brightness = 200), degrees = 90)
HE_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA7/HE_landmarks_corrected.csv", row.names = 1) 
Xenium_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA7/Xenium_landmarks_corrected.csv", row.names = 1) 
imageinfo_query <- magick::image_info(query_image)
imageinfo_ref <- magick::image_info(ref_image)

# transform validation markers (auto)
mapping_parameters <- readRDS("parameters/lungTMA7_HE_vs_DAPI_homography.rds")
Xenium_landmarks_auto <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_auto$y <- imageinfo_query$height - Xenium_landmarks_auto$y
Xenium_landmarks_auto[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_auto), mapping_parameters$mapping$`2`)
Xenium_landmarks_auto$y <- imageinfo_ref$height - Xenium_landmarks_auto$y
diff_auto <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_auto[,c("x", "y")])^2))

# transform validation markers (affine)
mapping_parameters <- readRDS("parameters/lungTMA7_HE_vs_DAPI_affine.rds")
Xenium_landmarks_affine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine$y <- imageinfo_query$height - Xenium_landmarks_affine$y
Xenium_landmarks_affine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine$y <- imageinfo_ref$height - Xenium_landmarks_affine$y
diff_affine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine[,c("x", "y")])^2))

# transform validation markers (affine + TPS)
mapping_parameters <- readRDS("parameters/lungTMA7_HE_vs_DAPI_affinetps.rds")
Xenium_landmarks_affinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affinetps$y <- imageinfo_query$height - Xenium_landmarks_affinetps$y
Xenium_landmarks_affinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_affinetps$y <- imageinfo_ref$height - Xenium_landmarks_affinetps$y
diff_affinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affinetps[,c("x", "y")])^2))

# transform validation markers (manual homo + tps)
mapping_parameters <- readRDS("parameters/lungTMA7_HE_vs_DAPI_manual_Affine.rds")
Xenium_landmarks_manualaffine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffine$y <- imageinfo_query$height - Xenium_landmarks_manualaffine$y
Xenium_landmarks_manualaffine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffine), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffine$y <- imageinfo_ref$height - Xenium_landmarks_manualaffine$y
diff_manualaffine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffine[,c("x", "y")])^2))

# transform validation markers (manual tps)
mapping_parameters <- readRDS("parameters/lungTMA7_HE_vs_DAPI_manual_TPS.rds")
Xenium_landmarks_manualtps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualtps$y <- imageinfo_query$height - Xenium_landmarks_manualtps$y
Xenium_landmarks_manualtps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualtps), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualtps$y <- imageinfo_ref$height - Xenium_landmarks_manualtps$y
diff_manualtps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualtps[,c("x", "y")])^2))

# transform validation markers (manual affine + tps)
mapping_parameters <- readRDS("parameters/lungTMA7_HE_vs_DAPI_manual_Affine_TPS.rds")
Xenium_landmarks_manualaffinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffinetps$y <- imageinfo_query$height - Xenium_landmarks_manualaffinetps$y
Xenium_landmarks_manualaffinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffinetps$y <- imageinfo_ref$height - Xenium_landmarks_manualaffinetps$y
diff_manualaffinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffinetps[,c("x", "y")])^2))

# transform validation markers (affine + BSpline)
mapping <- readRDS("parameters/lungTMA7_HE_vs_DAPI_affine.rds")
mapping$Method <- "Affine + Non-Rigid"
mapping$nonrigid <- "BSpline (SimpleITK)"
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(DAPI_image_list[[7]], brightness = 200)),
                               mapping_parameters = mapping, interactive = TRUE)
mapping_parameters <- xen_reg$mapping_parameters
Xenium_landmarks_affine_nr <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine_nr$y <- imageinfo_query$height - Xenium_landmarks_affine_nr$y
Xenium_landmarks_affine_nr[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine_nr), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine_nr$y <- imageinfo_ref$height - Xenium_landmarks_affine_nr$y
diff_affine_nr <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine_nr[,c("x", "y")])^2))

# compare markers
results <- data.frame(HE_landmarks[,c("x","y")], diff_manualaffine, diff_manualtps, diff_manualaffinetps, diff_affine, diff_affinetps, diff_auto, diff_affine_nr)
colnames(results) <- c('x_landmark', 'y_landmark', "diff_manualaffine", "diff_manualtps", "diff_manualaffinetps", "diff_affine", "diff_affinetps", "diff_auto", "diff_affine_nr")
dir.create("results/lungTMA7/")
write.csv(results, file = "results/lungTMA7/all_landmarks.csv", row.names = FALSE, quote = FALSE)

####
# Lung TMA 8 ####
####

####
## get parameters ####
####

ref_image <- HE_image_list[[8]]
query_image <- DAPI_image_list[[8]]
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA8_HE_vs_DAPI_homography.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA8_HE_vs_DAPI_homography.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA8_HE_vs_DAPI_affine.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA8_HE_vs_DAPI_affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(query_image, brightness = 200)), 
                               mapping_parameters = readRDS("parameters/lungTMA8_HE_vs_DAPI_affinetps.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA8_HE_vs_DAPI_affinetps.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA8_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA8_HE_vs_DAPI_manual_Affine.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA8_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA8_HE_vs_DAPI_manual_TPS.rds")

xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_rotate(image_modulate(query_image, brightness = 200), degree = 90)), 
                               mapping_parameters = readRDS("../../../../ImageAlignmentBenchmark/scripts/landmark_selection/landmarks/voltron/lungTMA8_affine_transformation_points.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "parameters/lungTMA8_HE_vs_DAPI_manual_Affine_TPS.rds")

####
## transformation of validation landmarks ####
####

# landmarks and images
ref_image <- HE_image_list[[8]]
query_image <- magick::image_rotate(image_modulate(DAPI_image_list[[8]], brightness = 200), degrees = 90)
HE_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA8/HE_landmarks_corrected.csv", row.names = 1) 
Xenium_landmarks <- read.csv("../../landmarks/LungTMA/lungTMA8/Xenium_landmarks_corrected.csv", row.names = 1) 
imageinfo_query <- magick::image_info(query_image)
imageinfo_ref <- magick::image_info(ref_image)

# transform validation markers (auto)
mapping_parameters <- readRDS("parameters/lungTMA8_HE_vs_DAPI_homography.rds")
Xenium_landmarks_auto <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_auto$y <- imageinfo_query$height - Xenium_landmarks_auto$y
Xenium_landmarks_auto[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_auto), mapping_parameters$mapping$`2`)
Xenium_landmarks_auto$y <- imageinfo_ref$height - Xenium_landmarks_auto$y
diff_auto <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_auto[,c("x", "y")])^2))

# transform validation markers (affine)
mapping_parameters <- readRDS("parameters/lungTMA8_HE_vs_DAPI_affine.rds")
Xenium_landmarks_affine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine$y <- imageinfo_query$height - Xenium_landmarks_affine$y
Xenium_landmarks_affine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine$y <- imageinfo_ref$height - Xenium_landmarks_affine$y
diff_affine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine[,c("x", "y")])^2))

# transform validation markers (affine + TPS)
mapping_parameters <- readRDS("parameters/lungTMA8_HE_vs_DAPI_affinetps.rds")
Xenium_landmarks_affinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affinetps$y <- imageinfo_query$height - Xenium_landmarks_affinetps$y
Xenium_landmarks_affinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_affinetps$y <- imageinfo_ref$height - Xenium_landmarks_affinetps$y
diff_affinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affinetps[,c("x", "y")])^2))

# transform validation markers (manual homo + tps)
mapping_parameters <- readRDS("parameters/lungTMA8_HE_vs_DAPI_manual_Affine.rds")
Xenium_landmarks_manualaffine <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffine$y <- imageinfo_query$height - Xenium_landmarks_manualaffine$y
Xenium_landmarks_manualaffine[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffine), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffine$y <- imageinfo_ref$height - Xenium_landmarks_manualaffine$y
diff_manualaffine <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffine[,c("x", "y")])^2))

# transform validation markers (manual tps)
mapping_parameters <- readRDS("parameters/lungTMA8_HE_vs_DAPI_manual_TPS.rds")
Xenium_landmarks_manualtps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualtps$y <- imageinfo_query$height - Xenium_landmarks_manualtps$y
Xenium_landmarks_manualtps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualtps), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualtps$y <- imageinfo_ref$height - Xenium_landmarks_manualtps$y
diff_manualtps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualtps[,c("x", "y")])^2))

# transform validation markers (manual affine + tps)
mapping_parameters <- readRDS("parameters/lungTMA8_HE_vs_DAPI_manual_Affine_TPS.rds")
Xenium_landmarks_manualaffinetps <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_manualaffinetps$y <- imageinfo_query$height - Xenium_landmarks_manualaffinetps$y
Xenium_landmarks_manualaffinetps[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_manualaffinetps), mapping_parameters$mapping$`2`)
Xenium_landmarks_manualaffinetps$y <- imageinfo_ref$height - Xenium_landmarks_manualaffinetps$y
diff_manualaffinetps <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_manualaffinetps[,c("x", "y")])^2))

# transform validation markers (affine + BSpline)
mapping <- readRDS("parameters/lungTMA8_HE_vs_DAPI_affine.rds")
mapping$Method <- "Affine + Non-Rigid"
mapping$nonrigid <- "BSpline (SimpleITK)"
xen_reg <- registerSpatialData(reference_spatdata = importImageData(ref_image),
                               query_spatdata = importImageData(image_modulate(DAPI_image_list[[8]], brightness = 200)),
                               mapping_parameters = mapping, interactive = TRUE)
mapping_parameters <- xen_reg$mapping_parameters
Xenium_landmarks_affine_nr <- Xenium_landmarks[,c("x","y")]
Xenium_landmarks_affine_nr$y <- imageinfo_query$height - Xenium_landmarks_affine_nr$y
Xenium_landmarks_affine_nr[,c("x", "y")] <- applyMapping(as.matrix(Xenium_landmarks_affine_nr), mapping_parameters$mapping$`2`)
Xenium_landmarks_affine_nr$y <- imageinfo_ref$height - Xenium_landmarks_affine_nr$y
diff_affine_nr <- sqrt(rowSums((HE_landmarks[,c("x","y")] - Xenium_landmarks_affine_nr[,c("x", "y")])^2))

# compare markers
results <- data.frame(HE_landmarks[,c("x","y")], diff_manualaffine, diff_manualtps, diff_manualaffinetps, diff_affine, diff_affinetps, diff_auto, diff_affine_nr)
colnames(results) <- c('x_landmark', 'y_landmark', "diff_manualaffine", "diff_manualtps", "diff_manualaffinetps", "diff_affine", "diff_affinetps", "diff_auto", "diff_affine_nr")
dir.create("results/lungTMA8/")
write.csv(results, file = "results/lungTMA8/all_landmarks.csv", row.names = FALSE, quote = FALSE)
