# library(VoltRon)
library(Seurat)
library(harmony)
library(dplyr)
library(ComplexHeatmap)
library(xlsx)

####
## Import Spatial Data ####
####

# import
vr2 <- importXenium("COVID19_lung_xenium/outs",
                    resolution_level = 3, overwrite_resolution = FALSE, import_molecules = TRUE)

# visualize molecules 
vrSpatialPlot(vr2, assay = "Xenium_mol", group.by = "gene", group.ids = c("S2_N", "S2_orf1ab"), n.tile = 900, alpha = 0.7)

####
## Subset all sections ####
####

# get subsets
subsetting_all <- subset(vr2, interactive = TRUE)

# use subset information and subset
vr2_list <- list()
samples <- c("prolonged case 4", "acute case 3", "control case 2", "acute case 1", "acute case 2",
             "prolonged case 5", "prolonged case 3", "control case 1")
for(i in 1:length(subsetting_all)){
  print(samples[i])
  vr2_list[[i]] <- subsetting_all[[i]]
  vr2_list[[i]]$Sample <- samples[i]
}
vr2_merged <- merge(vr2_list[[1]], vr2_list[-1])

####
## Analyze ####
####

vr2_merged_acute1 <- subset(vr2_merged, samples = "acute case 1")

####
### Processing and Embeddings ####
####

# processing
vr2_merged_acute1 <- subset(vr2_merged_acute1, Count > 5)
vr2_merged_acute1 <- normalizeData(vr2_merged_acute1, method = "LogNorm")

# embedding
selected_features <- vrFeatures(vr2_merged_acute1)
selected_features <- selected_features[!selected_features %in% c("S2_N", "S2_orf1ab")]
vr2_merged_acute1 <- getPCA(vr2_merged_acute1, features = selected_features, dims = 15, overwrite = TRUE)
vr2_merged_acute1 <- getUMAP(vr2_merged_acute1, dims = 1:15, overwrite = TRUE)

# visualization
vrEmbeddingPlot(vr2_merged_acute1, embedding = "umap")

####
### Clustering ####
####

# clusters
vr2_merged_acute1 <- getProfileNeighbors(vr2_merged_acute1, dims = 1:15, k = 10, method = "SNN", data.type = "pca", graph.key = "SNN")
for(i in c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)){
  vr2_merged_acute1 <- getClusters(vr2_merged_acute1, resolution = i, label = paste0("Cluster_", i), graph = "SNN")
}

# visualize
g_list <- list()
g_list[["Cluster_0.7"]] <- vrEmbeddingPlot(vr2_merged_acute1, embedding = "umap", group.by = "Cluster_0.7", label = T)
g_list[["Cluster_0.8"]] <- vrEmbeddingPlot(vr2_merged_acute1, embedding = "umap", group.by = "Cluster_0.8", label = T)
g_list[["Cluster_0.9"]] <- vrEmbeddingPlot(vr2_merged_acute1, embedding = "umap", group.by = "Cluster_0.9", label = T)
g_list[["Cluster_1"]] <- vrEmbeddingPlot(vr2_merged_acute1, embedding = "umap", group.by = "Cluster_1", label = T)
g_list[["Cluster_1.1"]] <- vrEmbeddingPlot(vr2_merged_acute1, embedding = "umap", group.by = "Cluster_1.1", label = T)
g_list[["Cluster_1.2"]] <- vrEmbeddingPlot(vr2_merged_acute1, embedding = "umap", group.by = "Cluster_1.2", label = T)
g_list[["Cluster_1.3"]] <- vrEmbeddingPlot(vr2_merged_acute1, embedding = "umap", group.by = "Cluster_1.3", label = T)
ggarrange(plotlist = g_list, ncol = 3, nrow = 3)

####
### Marker Analysis ####
####

####
#### res 1.3 ####
####

# save as Seurat
vr2_merged_acute1_seu <- VoltRon::as.Seurat(vr2_merged_acute1, cell.assay = "Xenium", type = "image")

# marker analysis
vr2_merged_acute1_seu <- NormalizeData(vr2_merged_acute1_seu, scale.factor = 100)
Idents(vr2_merged_acute1_seu) <- "Cluster_1.3"
markers <- FindAllMarkers(vr2_merged_acute1_seu)
topmarkers <- markers %>% 
  group_by(cluster) %>%
  filter(avg_log2FC > 0.5, pct.1 > 0.2, p_val_adj < 0.05) %>%
  slice_max(order_by = avg_log2FC, n = 10)

# add xenium annotations
xenium_annotations <- read.xlsx(file = "../data/xenium_annnotation.xlsx", sheetName = "Sheet1")
markers <- markers %>% left_join(xenium_annotations)
topmarkers <- topmarkers %>% left_join(xenium_annotations)
topmarkers$annotation.xenium[is.na(topmarkers$annotation.xenium)] <- "undefined"
topmarkers <- topmarkers %>%
  group_by(cluster) %>%
  mutate(FinalAnnotation = names(table(annotation.xenium))[which.max(table(annotation.xenium))])

# heatmap 
topmarkers <- markers %>% 
  group_by(cluster) %>%
  filter(avg_log2FC > 0.5, pct.1 > 0.1, p_val_adj < 0.05) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>% 
  arrange(match(cluster, levels(factor(vr2_merged_acute1$cluster))))

# heatmap
selected_features <- unique(c(topmarkers$gene,"S2_N", "S2_orf1ab"))
vrHeatmapPlot(vr2_merged_acute1, assay = "Xenium", features = selected_features, group.by = "Cluster_1.3", show_row_names = TRUE, cluster_rows = FALSE)

# violin plot
vrEmbeddingFeaturePlot(vr2_merged_acute1, assay = "Xenium", features = "COL4A1", embedding = "umap")
vrViolinPlot(vr2_merged_acute1, assay = "Xenium", features = "KCNK3", group.by = "Cluster_1.3")

####
### Annotation ####
####

topmarkers_celltype <- topmarkers %>%
  select(cluster, FinalAnnotation) %>%
  distinct()

# 1       Macrophages_1             
# 2       undefined (prev: Ciliated cells)             
# 3       Fibroblasts                
# 4       AT1&2 cells      
# 5       undefined 
# 6       Endothelial cells          
# 7       Macrophages_2                
# 8       Lymphatic endothelial cells
# 9       Classical monocytes        
# 10      Smooth muscle cells        
# 11      T cells                    
# 12      undefined    
# 13      NK cells  
# 14      Myofibroblasts                
# 15      Capillary cells_1            
# 16      Pericyte_2                  
# 17      undefined (prev: Capillary cells_2)            

# annotations
annotations <- c("Macrophages_1", "undefined", "Fibroblasts", "AT1&2 cells", "undefined", "Endothelial cells ", "Macrophages_2", "Lymphatic endothelial cells", 
                 "Classical monocytes", "Smooth muscle cells", "T cells", "undefined", "NK cells", "Myofibroblasts", "Capillary cells_1", "Pericyte_2", "undefined")
# annotate cells
vr2_merged_acute1$CellType <- annotations[match(vr2_merged_acute1$Cluster_1.3, 1:length(annotations))]
vrEmbeddingPlot(vr2_merged_acute1, embedding = "umap", group.by = "CellType", label = T)

####
### Subcluster undefined ####
####

vr2_merged_acute1_undefined <- subset(vr2_merged_acute1, subset = CellType == "undefined")

# processing
vr2_merged_acute1_undefined <- normalizeData(vr2_merged_acute1_undefined, method = "LogNorm")
selected_features <- vrFeatures(vr2_merged_acute1_undefined)
selected_features <- selected_features[!selected_features %in% c("S2_N", "S2_orf1ab")]
vr2_merged_acute1_undefined <- getPCA(vr2_merged_acute1_undefined, features = selected_features, dims = 10, type = "pca_undefined")
vr2_merged_acute1_undefined <- getUMAP(vr2_merged_acute1_undefined, dims = 1:10, data.type = "pca_undefined", umap.key = "umap_undefined")
vrEmbeddingPlot(vr2_merged_acute1_undefined, embedding = "umap_undefined")

# clusters q q
vr2_merged_acute1_undefined <- getProfileNeighbors(vr2_merged_acute1_undefined, dims = 1:10, k = 10, method = "SNN", data.type = "pca_undefined", graph.key = "SNN_undefined")
for(i in c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)){
  vr2_merged_acute1_undefined <- getClusters(vr2_merged_acute1_undefined, resolution = i, label = paste0("Cluster_", i), graph = "SNN_undefined")
}

# visualize
g_list <- list()
g_list[["Cluster_0.2"]] <- vrEmbeddingPlot(vr2_merged_acute1_undefined, embedding = "umap_undefined", group.by = "Cluster_0.2", label = T)
g_list[["Cluster_0.4"]] <- vrEmbeddingPlot(vr2_merged_acute1_undefined, embedding = "umap_undefined", group.by = "Cluster_0.4", label = T)
g_list[["Cluster_0.6"]] <- vrEmbeddingPlot(vr2_merged_acute1_undefined, embedding = "umap_undefined", group.by = "Cluster_0.6", label = T)
g_list[["Cluster_0.8"]] <- vrEmbeddingPlot(vr2_merged_acute1_undefined, embedding = "umap_undefined", group.by = "Cluster_0.8", label = T)
g_list[["Cluster_1"]] <- vrEmbeddingPlot(vr2_merged_acute1_undefined, embedding = "umap_undefined", group.by = "Cluster_1", label = T)
g_list[["Cluster_1.2"]] <- vrEmbeddingPlot(vr2_merged_acute1_undefined, embedding = "umap_undefined", group.by = "Cluster_1.2", label = T)
g_list[["Cluster_1.4"]] <- vrEmbeddingPlot(vr2_merged_acute1_undefined, embedding = "umap_undefined", group.by = "Cluster_1.4", label = T)
ggarrange(plotlist = g_list, ncol = 3, nrow = 3)

####
#### res 1.2 ####
####

# marker analysis
vr2_merged_acute1_undefined_seu <- VoltRon::as.Seurat(vr2_merged_acute1_undefined, cell.assay = "Xenium", type = "image")
vr2_merged_acute1_undefined_seu <- NormalizeData(vr2_merged_acute1_undefined_seu, scale.factor = 10000)
Idents(vr2_merged_acute1_undefined_seu) <- "Cluster_1.2"
markers <- FindAllMarkers(vr2_merged_acute1_undefined_seu)
topmarkers <- markers %>% 
  group_by(cluster) %>%
  filter(avg_log2FC > 0.5, pct.1 > 0.2, p_val_adj < 0.05) %>%
  slice_max(order_by = avg_log2FC, n = 10)
markers <- markers %>% left_join(xenium_annotations)
topmarkers <- topmarkers %>% left_join(xenium_annotations)
topmarkers$annotation.xenium[is.na(topmarkers$annotation.xenium)] <- "undefined"
topmarkers <- topmarkers %>%
  group_by(cluster) %>%
  mutate(FinalAnnotation = names(table(annotation.xenium))[which.max(table(annotation.xenium))])

# embedding
vrEmbeddingPlot(vr2_merged_acute1_undefined, embedding = "umap_undefined", group.by = "Cluster_1.2", label = T)

# check markers
vrEmbeddingFeaturePlot(vr2_merged_acute1_undefined, embedding = "umap_undefined", features =  c("S2_N", "S2_orf1ab"))
vrViolinPlot(vr2_merged_acute1_undefined, features =  c("S2_N", "S2_orf1ab"), group.by = "Cluster_0.8")
vrEmbeddingFeaturePlot(vr2_merged_acute1_undefined, embedding = "umap_undefined", features =  c("SCGB3A2", "MUC1", "EPCAM"))

vrEmbeddingFeaturePlot(vr2_merged_acute1_undefined, features =  c("TMC5", "SCGB3A2"), embedding = "umap_undefined")

####
#### Move annotations to main object ####
####

# annotation
vr2_merged_acute1_undefined$CellType <- paste0(vr2_merged_acute1_undefined$CellType, "_", vr2_merged_acute1_undefined$Cluster_1.2)

# move annotations to main object
temp <- vr2_merged_acute1$CellType
names(temp) <- vrSpatialPoints(vr2_merged_acute1)
temp[vrSpatialPoints(vr2_merged_acute1_undefined)] <- vr2_merged_acute1_undefined$CellType
vr2_merged_acute1$CellType <- temp

# check markers
vrViolinPlot(vr2_merged_acute1, features =  c("C1R", "C1S"), group.by = "CellType")
vrEmbeddingPlot(vr2_merged_acute1, embedding = "umap", group.by = "CellType", label = T)

# check markers
vr2_merged_acute1_seu <- VoltRon::as.Seurat(vr2_merged_acute1, cell.assay = "Xenium", type = "image")
vr2_merged_acute1_seu <- NormalizeData(vr2_merged_acute1_seu, scale.factor = 100)
Idents(vr2_merged_acute1_seu) <- "CellType"
markers <- FindAllMarkers(vr2_merged_acute1_seu)
topmarkers <- markers %>% 
  group_by(cluster) %>%
  filter(avg_log2FC > 0.5, pct.1 > 0.2, p_val_adj < 0.05) %>%
  slice_max(order_by = avg_log2FC, n = 10)
markers <- markers %>% left_join(xenium_annotations)
topmarkers <- topmarkers %>% left_join(xenium_annotations)
topmarkers$annotation.xenium[is.na(topmarkers$annotation.xenium)] <- "undefined"
topmarkers <- topmarkers %>%
  group_by(cluster) %>%
  mutate(FinalAnnotation = names(table(annotation.xenium))[which.max(table(annotation.xenium))])
topmarkers <- markers %>% 
  group_by(cluster) %>%
  filter(avg_log2FC > 0.5, pct.1 > 0.2, p_val_adj < 0.05) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>% 
  arrange(match(cluster, levels(factor(vr2_merged_acute1$cluster))))
selected_features <- unique(c(topmarkers$gene,"S2_N", "S2_orf1ab"))
vrHeatmapPlot(vr2_merged_acute1, assay = "Xenium", features = selected_features, group.by = "CellType", show_row_names = TRUE, cluster_rows = FALSE)

####
### Final Annotation ####
####

# Cluster_1.2
# 1 = S2++ cells
# 2 = undefined
# 3 = undefined
# 4 = S2++ cells
# 5 = Ciliated/Club Cells
# 6 = Myofibroblasts
# 7 = Myofibroblasts
# 8 = S2++ cells
# 9 = undefined

# annotation
vr2_merged_acute1$CellType[vr2_merged_acute1$CellType %in% c("undefined_1", "undefined_4", "undefined_8")] <- "S2 cells"
vr2_merged_acute1$CellType[vr2_merged_acute1$CellType %in% c("undefined_2", "undefined_3", "undefined_9")] <- "undefined"
vr2_merged_acute1$CellType[vr2_merged_acute1$CellType %in% c("undefined_6", "undefined_7")] <- "Myofibroblasts"
vr2_merged_acute1$CellType[vr2_merged_acute1$CellType %in% c("undefined_5")] <- "Ciliated/Club Cells"
vrEmbeddingPlot(vr2_merged_acute1, embedding = "umap", group.by = "CellType", label = T)

# new data
vr2_merged_acute1_comp <- subset(vr2_merged, samples = "acute case 1")
vr2_merged_acute1_comp$CellType <- "undefined"
vr2_merged_acute1_comp$CellType[match(vrSpatialPoints(vr2_merged_acute1), vrSpatialPoints(vr2_merged_acute1_comp))] <- vr2_merged_acute1$CellType
vrEmbeddings(vr2_merged_acute1_comp, type = "pca") <- vrEmbeddings(vr2_merged_acute1, type = "pca")
vrEmbeddings(vr2_merged_acute1_comp, type = "umap") <- vrEmbeddings(vr2_merged_acute1, type = "umap")
vrEmbeddingPlot(vr2_merged_acute1_comp, embedding = "umap", group.by = "CellType", label = T)
vrSpatialPlot(vr2_merged_acute1_comp, assay = "Xenium_mol", group.by = "gene", 
              group.ids = c("S2_N", "S2_orf1ab"), interactive = FALSE, pt.size = 0.1, n.tile = 300)

# save annotated data
vr2_merged_acute1_comp$CellType[vr2_merged_acute1_comp$CellType == "Pericyte_2"] <- "Pericytes"
vr2_merged_acute1_comp$CellType[vr2_merged_acute1_comp$CellType == "Capillary cells_1"] <- "Capillary cells"
vr2_merged_acute1_comp$CellType[vr2_merged_acute1_comp$CellType == "S2 cells"] <- "H.I. Cells"
saveRDS(vr2_merged_acute1_comp, file = "data/acutecase1_annotated.rds")

####
## Downstream Visualization ####
####

# read data
vr2_merged_acute1 <- readRDS(file = "data/acutecase1_annotated.rds")
vr2_merged_acute1_comp <- subset(vr2_merged, samples = "acute case 1")

vrMainAssay(vr2_merged_acute1) <- "Xenium"
vrSpatialPlot(vr2_merged_acute1, group.by = "CellType", group.ids = c("Infected cells"))
vrSpatialFeaturePlot(vr2_merged_acute1, features = c("S2_N", "S2_orf1ab"), alpha = 1)
vrViolinPlot(vr2_merged_acute1, group.by = "CellType", features = c("S2_N", "S2_orf1ab"))

vrSpatialPlot(vr2_merged_acute1, assay = "Xenium_mol", group.by = "gene", group.ids = c("S2_N", "S2_orf1ab"), n.tile = 300, 
              background = c("main", "H&E"), alpha = 0.7)

####
## Register images ####
####

# get image
imgdata <- importImageData("../../../data/HelenaAnja_13092023/HE_stainings/middlerowmiddle.jpg", 
                           tile.size = 10, 
                           segments = "../../../data/HelenaAnja_13092023/HE_stainings/annotations/acutecase1_membrane.geojson", 
                           sample_name = "acute case 1 (HE)")
imgdata <- flipCoordinates(imgdata, assay = "ROIAnnotation")
imgdata

# visualization
vrSpatialPlot(imgdata, assay = "ROIAnnotation", group.by = "Sample", alpha = 0.7, interactive = FALSE)

# get subset
vr2_merged_acute1 <- modulateImage(vr2_merged_acute1, brightness = 300, channel = "DAPI")

# width 1859
xen_reg <- registerSpatialData(object_list = list(vr2_merged_acute1, imgdata))
imgdata_reg <- xen_reg$registered_spat[[2]]
vrSpatialPlot(imgdata_reg, assay = "ROIAnnotation", group.by = "Sample", alpha = 0.5, background = c("image_1_reg"))

# add assays and images
vrMainAssay(imgdata_reg) <- "ROIAnnotation"
vrImages(vr2_merged_acute1[["Assay7"]], name = "main", channel = "H&E") <- vrImages(imgdata_reg, assay = "Assay1", name = "image_1_reg")
vrImages(vr2_merged_acute1[["Assay8"]], name = "main", channel = "H&E") <- vrImages(imgdata_reg, assay = "Assay1", name = "image_1_reg")
vr2_merged_acute1 <- addAssay(vr2_merged_acute1,
                       assay = imgdata_reg[["Assay2"]],
                       metadata = Metadata(imgdata_reg, assay = "ROIAnnotation"),
                       assay_name = "ROIAnnotation",
                       sample = "acute case 1", layer = "Section1")

# visualize 
vrSpatialPlot(vr2_merged_acute1, assay = "Assay7", plot.segments = TRUE, background = c("main","H&E"), alpha = 0.5)
vrSpatialPlot(vr2_merged_acute1, assay = "Assay17", plot.segments = TRUE, background = "image_1_reg", alpha = 0.5)

# transfer data
vrMainAssay(vr2_merged_acute1) <- "ROIAnnotation"
vr2_merged_acute1$Region <- vrSpatialPoints(vr2_merged_acute1)
vrMainImage(vr2_merged_acute1[["Assay9"]]) <- "image_1_reg"
vr2_merged_acute1 <- transferData(object = vr2_merged_acute1, from = "Assay9", to = "Assay8", features = "Region")

# visualize transcripts
vrSpatialPlot(vr2_merged_acute1, assay = "Assay8", group.by = "Region", background = c("main","H&E"), alpha = 0.5, n.tile = 300)

####
## as.AnnData ####
####

# acute case 1
vr2_subset <- subset(vr2_merged_acute1, assays = "Assay7")
as.AnnData(vr2_subset, file = "data/acutecase1_annotated.h5ad", flip_coordinates = TRUE, channel = "H&E")
