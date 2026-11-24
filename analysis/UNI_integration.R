library(Seurat)
library(xlsx)
library(magick)

####
# Xenium Breast H&E ####
####

# import
Xen_R1 <- importXenium("../../../../../data/10X_Xenium_Visium/Xenium_R1/outs",
                       sample_name = "XeniumR1", resolution_level = 3,
                       overwrite_resolution = TRUE)
Xen_R1_image <- importImageData("../../../../../data/10X_Xenium_Visium/Xenium_R1/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image.tif",
                                sample_name = "XeniumR1",
                                channels = "H&E", tile.size = 128)

xenium_annotation <- read.csv("../../../../../test/mainWorkflows/data/BreastTile/UNI/Xenium_hBreast_v1_metadata.csv")

####
## Clustering ####
####

# get embed data
embed_data <- readRDS("../../../../../test/mainWorkflows/data/BreastTile/UNI/uni2h_features_processed.rds")
ind <- sapply(vrSpatialPoints(Xen_R1_image), function(x) strsplit(x, split = "_")[[1]][1])
embed_data <- embed_data[ind,]
rownames(embed_data) <- vrSpatialPoints(Xen_R1_image)

# add embedding
vrEmbeddings(Xen_R1_image, type = "dinov2", overwrite = TRUE) <-  embed_data
Xen_R1_image <- getPCA(Xen_R1_image, dims = 50, data.type = "dinov2")
Xen_R1_image <- getUMAP(Xen_R1_image)
Xen_R1_image <- getClusters(Xen_R1_image, method = "kmeans", label = "Cluster_kmeans", data.type = "pca", nclus = 13)

####
## Registration ####
####

# Register
# xen_reg <- registerSpatialData(object_list = list(Xen_R1_image, Xen_R1))
xen_reg <- registerSpatialData(object_list = list(Xen_R1_image, Xen_R1), 
                               mapping_parameters = readRDS("auxiliary/he_registration_parameters.rds"), 
                               interactive = FALSE)

####
## Data Transfer ####
####

# add assay
Xen_R1_reg <- xen_reg$registered_spat[[2]]
Xen_R1_merged <- Xen_R1_image
Xen_R1_merged <- addAssay(Xen_R1_merged, 
                          metadata = Metadata(Xen_R1_reg), 
                          assay = Xen_R1_reg[["Assay1"]],
                          assay_name = "Xenium",
                          sample = "XeniumR1")

# transfer
Xen_R1_merged <- transferData(Xen_R1_merged, 
                              from = "Assay1", 
                              to = "Assay2", 
                              features = "Cluster_kmeans", 
                              expand = FALSE)

# transfer image
vrImages(Xen_R1_merged[["Assay2"]], name = "main_reg", channel = "H&E") <- 
  vrImages(Xen_R1_merged, assay = "Assay1")

vrSpatialPlot(Xen_R1_merged, group.by = "Cluster_kmeans")
vrSpatialPlot(Xen_R1_merged, group.by = "Cluster_kmeans", assay = "Xenium")

####
## Marker analysis ####
####

Xen_data_seu <- VoltRon::as.Seurat(Xen_R1_merged, cell.assay = "Xenium", type = "image")
Xen_data_seu <- NormalizeData(Xen_data_seu, scale.factor = 100)
Idents(Xen_data_seu) <- "Cluster_kmeans"
markers_image <- FindAllMarkers(Xen_data_seu)
markers_image <- markers_image %>% left_join(xenium_annotation, by = c("gene" = "Genes"))

# top markers
markers_image_top <- markers_image %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05, pct.1 > 0.5)

vrHeatmapPlot(Xen_R1_merged, group.by = "Cluster_kmeans", assay = "Xenium",
              features = unique(markers_image_top$gene), show_row_names = TRUE)

####
## Annotation ####
####

vrMainAssay(Xen_R1_merged) <- "Xenium"
annotation <- Xen_R1_merged$Cluster_kmeans 
annotation[annotation %in% c(11,13)] <- "Tumor"
annotation[annotation %in% c(8)] <- "Myoepithelial"
annotation[annotation %in% c(6)] <- "Immune"
annotation[annotation %in% c(7,10)] <- "Stroma"
annotation[annotation %in% c(1,2,3,4,5,9,12)] <- "Background"
Xen_R1_merged$annotation <- annotation

vrMainAssay(Xen_R1_merged) <- "ImageData"
annotation <- Xen_R1_merged$Cluster_kmeans 
annotation[annotation %in% c(11,13)] <- "Tumor"
annotation[annotation %in% c(8)] <- "Myoepithelial"
annotation[annotation %in% c(6)] <- "Immune"
annotation[annotation %in% c(7,10)] <- "Stroma"
annotation[annotation %in% c(1,2,3,4,5,9,12)] <- "Background"
Xen_R1_merged$annotation <- annotation

####
### Marker Analysis ####
####

Xen_data_seu <- VoltRon::as.Seurat(Xen_R1_merged, cell.assay = "Xenium", type = "image")
Xen_data_seu <- NormalizeData(Xen_data_seu, scale.factor = 100)
Idents(Xen_data_seu) <- "annotation"
markers_image <- FindAllMarkers(Xen_data_seu)
markers_image <- markers_image %>% left_join(xenium_annotation, by = c("gene" = "Genes"))

# top markers
markers_image_top <- markers_image %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05, pct.1 > 0.5)

vrHeatmapPlot(Xen_R1_merged, group.by = "annotation", assay = "Xenium",
              features = unique(markers_image_top$gene), show_row_names = TRUE)
vrSpatialPlot(Xen_R1_merged, group.by = "annotation")

####
## Investigate Markers ####
####

vrSpatialFeaturePlot(Xen_R1_merged, features = c("FASN", "CEACAM6"))

####
## Visualize ####
####

vrMainAssay(Xen_R1_merged) <- "Xenium"
Xen_R1_merged <- normalizeData(Xen_R1_merged, sizefactor = 1000)
Xen_R1_merged_subset <- subset(Xen_R1_merged, 
                               subset = annotation %in% unique(Xen_R1_merged$annotation)[unique(Xen_R1_merged$annotation) != "Background"])

# slide
vrSpatialPlot(Xen_R1_merged, assay = "ImageData", group.by = "annotation", alpha = 0.7, 
              colors = list(Background = "lightblue", Stroma = "brown1", Immune = "yellow", Tumor = "purple", Myoepithelial = "black")) + 
  labs(title = "") + 
  theme(legend.position = "none")

vrMainAssay(Xen_R1_merged) <- "ImageData"
Xen_R1_merged$annotation_new <- ifelse(Xen_R1_merged$annotation == "Stroma", 
                                       paste0("Stroma", Xen_R1_merged$Cluster_kmeans, sep = "_"),
                                       Xen_R1_merged$annotation)
vrSpatialPlot(Xen_R1_merged, assay = "ImageData", group.by = "annotation_new", alpha = 0.7) + 
  labs(title = "") + 
  theme(legend.position = "none")

# markers
vrMainAssay(Xen_R1_merged_subset) <- "Xenium"
vrHeatmapPlot(Xen_R1_merged_subset, group.by = "annotation", 
              features = unique(markers_image_top$gene), show_row_names = TRUE)
