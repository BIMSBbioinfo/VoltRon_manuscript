library(RBioFormats)
library(peakRAM)
library(bench)
library(magick)
library(VoltRon)

####
# Xenium breast 5000 ####
####

# import
Xen_breast <- importXenium("../../../../data/10X_Xenium_Visium/Xenium_R1/outs/", sample_name = "XeniumR1",
                       resolution_level = 2, overwrite_resolution = TRUE, import_molecules = TRUE)
ome.tiff <- "../../../../data/10X_Xenium_Visium/Xenium_R1/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image.tif"
read.metadata(ome.tiff)
Xen_breast_image <- importImageData(ome.tiff,
                                  sample_name = "XeniumImage",
                                  channel_names = "H&E", 
                                  tile.size = 100)

# register
# xen_reg <- registerSpatialData(object_list = c(Xen_breast_image, Xen_breast))
# saveRDS(xen_reg$mapping_parameters, file = "data/xenium_breast_parameters.rds")
mapping_parameters <- readRDS("data/xenium_breast_parameters.rds")
img <- vrImages(Xen_breast_image)
img2 <- vrImages(Xen_breast, assay = "Assay1")
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
  img_new <- VoltRon:::getRcppWarpImage(img, img2, mapping_parameters$mapping$`2`)
})

# mapping cell
coords <- vrCoordinates(Xen_breast)[,1:2]
res3 <- peakRAM({
  coords_new <- VoltRon:::applyMapping(coords, mapping_parameters$mapping$`2`)
})

# mapping molecules
coords_mol <- vrCoordinates(Xen_breast, assay = "Assay2")[,1:2]
res4 <- peakRAM({
  coords_mol_new <- VoltRon:::applyMapping(coords_mol, mapping_parameters$mapping$`2`)
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
  as.character(nrow(coords_mol)),
  as.character(res4$Elapsed_Time_sec),
  as.character(as_bench_bytes(res4$Peak_RAM_Used_MiB)),
  as.character(as_bench_bytes(object.size(image_data(img)) + object.size(image_data(img2)) + object.size(image_data(img_new))))
)
names(result) <- c("size_ref", "size_query", "time_reg", "mem_reg", "time_warp", "mem_warp", "n_cells", "time_map", "mem_map", 
                   "n_molecules", "time_map_mol", "mem_map_mol", "image_size")
write.csv(as.data.frame(result), file = "results/Xenium_breast.csv")
