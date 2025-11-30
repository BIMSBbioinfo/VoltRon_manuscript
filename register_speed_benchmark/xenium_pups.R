library(RBioFormats)
library(peakRAM)
library(bench)
library(magick)
library(VoltRon)

####
# Xenium pups ####
####

# import
Xen_pups <- importXenium("../../../../data/xenium/Xenium_Prime_Mouse_Pup_FFPE_outs/", sample_name = "XeniumR1",
                       resolution_level = 2, overwrite_resolution = FALSE, import_molecules = FALSE)
ome.tiff <- "../../../../data/xenium/Xenium_Prime_Mouse_Pup_FFPE_outs/Xenium_Prime_Mouse_Pup_FFPE_he_image.ome.tif"
read.metadata(ome.tiff)
Xen_pups_image <- importImageData("../../../../data/xenium/Xenium_Prime_Mouse_Pup_FFPE_outs/Xenium_Prime_Mouse_Pup_FFPE_he_image.ome.tif",
                                  sample_name = "XeniumImage", 
                                  tile.size = 100,
                                  resolution = 3, 
                                  series = 1)

# register
# xen_reg <- registerSpatialData(object_list = c(Xen_pups_image, Xen_pups),
#                                mapping_parameters = readRDS("data/xenium_pups_parameters.rds"))
# saveRDS(xen_reg$mapping_parameters, file = "data/xenium_pups_parameters.rds")
mapping_parameters <- readRDS("data/xenium_pups_parameters.rds")
img <- vrImages(Xen_pups_image)
img2 <- vrImages(Xen_pups)
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

# mapping
coords <- vrCoordinates(Xen_pups)[,1:2]
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
names(result) <- c("size_ref", "size_query", "time_reg", "mem_reg", "time_warp", "mem_warp", "n_cells", "time_map", "mem_map", "image_size")
write.csv(as.data.frame(result), file = "results/Xenium_pups.csv")
