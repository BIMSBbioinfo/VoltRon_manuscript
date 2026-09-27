library(RBioFormats)
library(peakRAM)
library(bench)
library(magick)
library(VoltRon)

####
# Xenium lung ####
####

# import Xenium
vr2 <- importXenium("../../../../data/HelenaAnja_13092023/output-XETG00046__0010726__Region_1__20230908__130559/",
                    resolution_level = 3, overwrite_resolution = TRUE, import_molecules = TRUE)
subset_info_list <- readRDS("../../../../analysis/HelenaAnja_13092023/data/subset_info_list.rds")[[1]]
Xen_lung <- subset(vr2, image = subset_info_list[4])

# import HE
HE_image <- magick::image_read("../../../../analysis/ImageAlignmentBenchmark/data/LungTMA/HE/lungTMA_4.jpg")
Xen_lung_image <- importImageData(HE_image,
                                  sample_name = "XeniumImage")

# register
# xen_reg <- registerSpatialData(object_list = c(Xen_lung_image, Xen_lung))
# saveRDS(xen_reg$mapping_parameters, file = "data/xenium_lung_parameters.rds")
mapping_parameters <- readRDS("data/xenium_lung_parameters.rds")
img <- vrImages(Xen_lung_image)
img2 <- vrImages(Xen_lung, assay = "Assay1")
img_list <- list(list(`H&E` = img), 
                 list(DAPI = img2))
start <- proc.time()
res1 <- peakRAM({
  xen_reg <- VoltRon:::computeAutomatedPairwiseTransform(image_list = img_list, 
                                               input = mapping_parameters, 
                                               channel_names = "art", ref = 1, 
                                               query = 2)
})

# image warping
res2 <- peakRAM({
  img_new <- VoltRon:::warpImage(img, img2, mapping_parameters$mapping$`2`)
})

# mapping cell
coords <- vrCoordinates(Xen_lung)[,1:2]
res3 <- peakRAM({
  coords_new <- VoltRon:::applyMapping(coords, mapping_parameters$mapping$`2`)
})

# mapping molecules
coords_mol <- vrCoordinates(Xen_lung, assay = "Assay2")[,1:2]
res4 <- peakRAM({
  coords_mol_new <- VoltRon:::applyMapping(coords_mol, mapping_parameters$mapping$`2`)
})

# collect results
result <- c(
  paste(VoltRon:::getImageInfo(img)[c("width", "height")], collapse = "x"),
  paste(VoltRon:::getImageInfo(img2)[c("width", "height")], collapse = "x"),
  as.character(res1$Elapsed_Time_sec),
  as.character(res1$Peak_RAM_Used_MiB),
  as.character(res2$Elapsed_Time_sec),
  as.character(res2$Peak_RAM_Used_MiB),
  as.character(nrow(coords)),
  as.character(res3$Elapsed_Time_sec),
  as.character(res3$Peak_RAM_Used_MiB),
  as.character(nrow(coords_mol)),
  as.character(res4$Elapsed_Time_sec),
  as.character(res4$Peak_RAM_Used_MiB),
  as.character(as_bench_bytes(object.size(image_data(img)) + object.size(image_data(img2)) + object.size(image_data(img_new))))
)
names(result) <- c("size_ref", "size_query", "time_reg", "mem_reg", "time_warp", "mem_warp", "n_cells", "time_map", "mem_map", 
                   "n_molecules", "time_map_mol", "mem_map_mol", "image_size")
write.csv(as.data.frame(result), file = "results/Xenium_lung.csv")
