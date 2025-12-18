library(RBioFormats)
library(peakRAM)
library(bench)
library(magick)
library(VoltRon)

####
# visium dlfpc ####
####

DLPFC_1 <- importVisium("../../../../data/10X_Visium_DLPFC/151673", sample_name = "DLPFC_1", resolution_level = "hires")
DLPFC_2 <- importVisium("../../../../data/10X_Visium_DLPFC/151674", sample_name = "DLPFC_2", resolution_level = "hires")
# DLPFC_1_subset <- subset(DLPFC_1, interactive = TRUE)
DLPFC_1_subset <- subset(DLPFC_1, image = "1323x1393+405+361")
# DLPFC_2_subset <- subset(DLPFC_2, interactive = TRUE)
DLPFC_2_subset <- subset(DLPFC_2, image = "1323x1393+405+361")

# register 1000 20
# xen_reg <- registerSpatialData(object_list = c(DLPFC_1_subset, DLPFC_2_subset))
# saveRDS(xen_reg$mapping_parameters, file = "data/visium_dlfpc_hires_parameters.rds")
mapping_parameters <- readRDS("data/visium_dlfpc_hires_parameters.rds")
img <- vrImages(DLPFC_1)
img2 <- vrImages(DLPFC_2)
img_list <- list(list(`H&E` = img), 
                 list(`H&E` = img2))
start <- proc.time()
res1 <- peakRAM({
  xen_reg <- VoltRon:::computeAutomatedPairwiseTransform(image_list = img_list, 
                                               input = mapping_parameters, 
                                               channel_names = "art", ref = 1, 
                                               query = 2)
})

# image warping
res2 <- peakRAM({
  img_new <- VoltRon:::getRcppWarpImage(img, img2, mapping_parameters$mapping$`2`)
})

# mapping
coords <- vrCoordinates(DLPFC_2)[,1:2]
res3 <- peakRAM({
  coords_new <- VoltRon:::applyMapping(coords, mapping_parameters$mapping$`2`)
})

# collect results
result <- c(
  paste(VoltRon:::getImageInfo(img)[c("width", "height")], collapse = "x"),
  paste(VoltRon:::getImageInfo(img2)[c("width", "height")], collapse = "x"),
  as.character(res1$Elapsed_Time_sec),
  as.character(as_bench_bytes(res1$Peak_RAM_Used_MiB)),
  as.character(res2$Elapsed_Time_sec),
  as.character(as_bench_bytes(res2$Peak_RAM_Used_MiB)),
  as.character(nrow(coords)),
  as.character(res3$Elapsed_Time_sec),
  as.character(as_bench_bytes(res3$Peak_RAM_Used_MiB)),
  as.character(as_bench_bytes(object.size(image_data(img)) + object.size(image_data(img2)) + object.size(image_data(img_new))))
)
names(result) <- c("size_ref", "size_query", "time_reg", "mem_reg", "time_warp", "mem_warp", "n_cells", "time_map", "mem_map", 
                   "image_size")
write.csv(as.data.frame(result), file = "results/visium_dlfpc_hires.csv")
