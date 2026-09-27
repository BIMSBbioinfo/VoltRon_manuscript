# library(VoltRon)
library(patchwork)
library(ggpubr)
library(mclust)
library(Seurat)
library(dplyr)
library(RBioFormats)
library(xlsx)
library(tiff)

####
# Import Xenium ####
####

# xenium 
# xen <- importXenium("Tonsil_xenium/outs", 
#                     resolution_level = 3, overwrite_resolution = TRUE, 
#                     sample_name = "XeniumR1")
xen <- importXenium("../../../../data/Erik/Erik_Alignment_TMA/out/output-XETG00420__0041844__Region_1__20250205__140630/",
                    resolution_level = 3, overwrite_resolution = TRUE,
                    sample_name = "XeniumR1")

# xenium annnotation
xenium_lung <- read.csv("../data/Tonsil/Xenium_hLung_v1_metadata.csv")

####
# Xenium Tonsil (Same Section) ####
####

####
## Import IF Data ####
####

# same section IF
# ome.tiff <- "Tonsil/Same section/Core_13.ome.tif"
ome.tiff <- "../../../../data/Erik/Erik_Alignment_TMA/out/Same section/Core_13.ome.tif"
ome_tiff_same_vr <- importImageData(ome.tiff, 
                                    tile.size = 100, 
                                    series = 1, 
                                    resolution = 1, 
                                    channels = c("DAPI", "cy5_cd20", "Cy5_CD21"), 
                                    is.RGB = FALSE)

####
## same section alignment ####
####

# Subset Xenium data
# xen <- subset(xen, interactive = TRUE)
# xen_subset <- xen$subsets[[1]]
xen_subset <- subset(xen, image = "5149x4945+7145+4380")

# register
xen_subset2 <- modulateImage(xen_subset, brightness = 500)
ome_tiff_same_vr2 <- modulateImage(ome_tiff_same_vr, brightness = 700, channel = "DAPI")
xen_reg <- registerSpatialData(object_list = list(ome_tiff_same_vr2, xen_subset2),
                               mapping_parameters = readRDS("../data/Tonsil/samesection_registration_parameters.rds"), 
                               interactive = TRUE)
xenium_reg <- xen_reg$registered_spat[[2]]

xen_reg <- registerSpatialData(object_list = list(ome_tiff_same_vr, xen_subset),
                               mapping_parameters = readRDS("../data/Tonsil/samesection_registration_parameters.rds"), 
                               interactive = TRUE)
# saveRDS(xen_reg$accuracy_metrics, file = "results/same_section_accuracy.rds")

# visualize accuracy
accuracy_metrics <- readRDS(file = "results/same_section_accuracy.rds")
VoltRon:::.plot_ssim_map(accuracy_metrics$ssim_maps$`2`$Coarse)

# visualize segments
vrImages(xenium_reg[["Assay1"]], name = "main_reg", channel = "CD20") <- vrImages(ome_tiff_same_vr, channel = "cy5_cd20")
vrImages(xenium_reg[["Assay1"]], name = "main_reg", channel = "CD21") <- vrImages(ome_tiff_same_vr, channel = "Cy5_CD21")

####
## Compare DAPI signals from both images ####
####

img1 <- vrImages(ome_tiff_same_vr, channel = "DAPI")
img2 <- vrImages(xenium_reg[["Assay1"]], name = "main_reg", channel = "DAPI")
img1_data <- magick::image_data(img1)
img1_data <- as.integer(img1_data)
img1_data <- array(img1_data, dim = dim(img1_data)[1:2])
img1_data <- as.matrix(img1_data)
img1_data <- img1_data/max(img1_data)
writeTIFF(img1_data, where = "results/img1.tiff", compression = "LZW")
img2_data <- magick::image_data(img2)
img2_data <- as.integer(img2_data)
img2_data <- array(img2_data, dim = dim(img2_data)[1:2])
img2_data <- as.matrix(img2_data)
img2_data <- img2_data/max(img2_data)
writeTIFF(img2_data, where = "results/img2.tiff", compression = "LZW")

img1 <- magick::image_modulate(img1, brightness = 700)
img2 <- magick::image_convert(img2, colorspace = "gray")
img1_data <- magick::image_data(img1)
img2_data <- magick::image_data(img2)
img1_data <- as.integer(img1_data)
img2_data <- as.integer(img2_data)
set.seed(1)
datax <- data.frame(img1 = as.vector(img1_data), img2 = as.vector(img2_data))
datax2 <- datax[sample(1:nrow(datax), 30000),]
ggplot(datax2, aes(x = img1, y = img2)) + 
  geom_point()
# plot(img1_data, img2_data)

####
## create and add aligned IF features ####
####

# get registered segments
segts <- vrSegments(flipCoordinates(xenium_reg))
generateGeoJSON(segts, "../data/Tonsil/polygons.geojson")

# add feature from QuPath
datax <- read.table("../data/Tonsil/measurements.txt", header = TRUE, sep = "\t")
rownames(datax) <- datax$Name
datax <- datax[vrSpatialPoints(xenium_reg),]
datax <- datax[,c(11,13)]
datax <- apply(datax, 2, function(x){
  x[is.na(x)] <- min(x, na.rm = TRUE)
  x
})
datax[,1] <- ifelse(datax[,1] < 1054, datax[,1],1054)
datax[,1] <- ifelse(datax[,1] > 331, datax[,1],331)
datax[,1] <- datax[,1] - min(datax[,1])
datax[,2] <- ifelse(datax[,2] < 1421, datax[,2],1421)
datax[,2] <- ifelse(datax[,2] > 491, datax[,2],491)
datax[,2] <- datax[,2] - min(datax[,2])
colnames(datax) <- c("CD20", "CD21")

# add to feature
datax <- t(datax)
xenium_reg <- addFeature(xenium_reg, data = datax, feature_name = "Protein")

####
## process and normalize ####
####

# subset
spatialpoints <- vrSpatialPoints(xenium_reg)[as.vector(Metadata(xenium_reg)$Count > 30)]
xenium_reg <- subset(xenium_reg, spatialpoints = spatialpoints)

# normalize RNA
xenium_reg <- normalizeData(xenium_reg, sizefactor = 1000)

# normalize Protein markers
vrMainFeatureType(xenium_reg) <- "Protein"
xenium_reg <- normalizeData(xenium_reg, method = "hyper.arcsine", scale = 0.2)

# visualize
vrMainFeatureType(xenium_reg) <- "Protein"
g1 <- vrSpatialFeaturePlot(xenium_reg, features = c("CD20", "CD21"), background.color = "black", n.tile = 300, ncol = 2, norm = TRUE)
vrMainFeatureType(xenium_reg) <- "RNA"
g2 <- vrSpatialFeaturePlot(xenium_reg, features = c("MS4A1", "CD38"), background.color = "black", n.tile = 300, ncol = 2)
g1 / g2

####
## on disk ####
####

xenium_reg <- saveVoltRon(xenium_reg, 
                          format = "HDF5VoltRon", 
                          output = "../data/Tonsil/Xenium_IF", 
                          replace = TRUE)
xenium_reg <- loadVoltRon(dir = "../data/Tonsil/Xenium_IF/")

####
## clustering (Xenium + IF) ####
####

# process and cluster
vrMainFeatureType(xenium_reg) <- "RNA"
xenium_reg <- getPCA(xenium_reg, 
                     feat_type = vrFeatureTypeNames(xenium_reg),
                     dims = 30)
xenium_reg <- getUMAP(xenium_reg, dims = 1:30)

# visualize features on embeddings 
vrMainFeatureType(xenium_reg) <- "Protein"
g1 <- vrEmbeddingFeaturePlot(xenium_reg, 
                             features = c("CD20", "CD21"), embedding = "umap", 
                             pt.size = 0.4)
vrMainFeatureType(xenium_reg) <- "RNA"
g2 <- vrEmbeddingFeaturePlot(xenium_reg, 
                             features = c("MS4A1", "CD38"), embedding = "umap", 
                             pt.size = 0.4)
g1 / g2

# clustering
xenium_reg <- getProfileNeighbors(xenium_reg, 
                                  dims = 1:30, method = "SNN")
for(res in c(0.7,0.9, 1.0, 1.1, 1.2, 1.3, 1.4)){
  print(res)
  xenium_reg <- getClusters(xenium_reg, resolution = res, 
                            label = paste("clusters", res, sep = "_"), graph = "SNN") 
}

# visualize clusters
gg_list <- list()
for(res in c(0.7,0.9, 1.0, 1.1, 1.2, 1.3, 1.4)){
  gg_list[[as.character(res)]] <-
    vrEmbeddingPlot(xenium_reg, 
                    group.by = paste("clusters", res, sep = "_"), 
                    embedding = "umap", 
                    pt.size = 0.4, label = TRUE)
}
ggpubr::ggarrange(plotlist = gg_list, ncol = 4, nrow = 2)

# visualize resolution 1
xenium_reg$Clusters <- xenium_reg$clusters_1
vrEmbeddingPlot(xenium_reg, 
                group.by = "Clusters",
                embedding = "umap", 
                pt.size = 0.4, label = TRUE)

# visualize cd21 expression across clusters
vrMainFeatureType(xenium_reg) <- "Protein"
vrViolinPlot(xenium_reg, features = "CD21", group.by = "Clusters")

# visualize cluster 14
vrSpatialPlot(xenium_reg, group.by = "Clusters", group.ids = 14,
              plot.segments = TRUE)

####
## Subclustering Germinal Center ####
####

spatialpoints <- vrSpatialPoints(xenium_reg)[xenium_reg$Clusters == 14]
xenium_reg_subset <- subset(xenium_reg, spatialpoints = spatialpoints)
vrSpatialPlot(xenium_reg_subset, group.by = "Clusters")
vrEmbeddingPlot(xenium_reg_subset, group.by = "Clusters", embedding = "umap")
xenium_reg_subset <- getClusters(xenium_reg_subset, resolution = 0.3, 
                                 label = "SubClusters", graph = "SNN") 
tmp <- xenium_reg_subset$SubClusters
tmp[tmp %in% c(1,2)] <- 1
tmp[tmp %in% c(3,4,5)] <- 2
xenium_reg_subset$SubClusters <- tmp
vrEmbeddingPlot(xenium_reg_subset, group.by = "SubClusters", embedding = "umap")
vrSpatialPlot(xenium_reg_subset, group.by = "SubClusters", plot.segments = TRUE)
vrMainFeatureType(xenium_reg_subset) <- "Protein"
vrEmbeddingFeaturePlot(xenium_reg_subset, features = "CD21", embedding = "umap")
vrMainFeatureType(xenium_reg_subset) <- "RNA"
vrEmbeddingFeaturePlot(xenium_reg_subset, features = "MKI67", embedding = "umap")

tmp <- setNames(xenium_reg$Clusters, vrSpatialPoints(xenium_reg))
tmp[vrSpatialPoints(xenium_reg_subset)] <- paste0(xenium_reg_subset$Clusters, "_", xenium_reg_subset$SubClusters)
names(tmp) <- NULL
xenium_reg$Clusters_withsub <- tmp

####
## marker analysis ####
####

# marker analysis
datax <- vrData(xenium_reg, feat_type = vrFeatureTypeNames(xenium_reg), norm = TRUE)
xenium_reg_seu <- CreateSeuratObject(as(datax, "dgCMatrix"), meta.data = as.data.frame(Metadata(xenium_reg)))
xenium_reg_seu$clusters_1 <- as.character(xenium_reg_seu$clusters_1)
new_datax <- datax
rownames(new_datax) <- rownames(xenium_reg_seu)
xenium_reg_seu <- SetAssayData(
  object = xenium_reg_seu,
  layer = "data",
  new.data = as(new_datax,"dgCMatrix"),
  assay = "RNA"
)

# xenium_reg_seu <- NormalizeData(xenium_reg_seu, scale.factor = 1000)
# Idents(xenium_reg_seu) <- "clusters_1"
Idents(xenium_reg_seu) <- "Clusters_withsub"
markers <- FindAllMarkers(xenium_reg_seu, features = Features(xenium_reg_seu))
markers$gene_symbol <- sapply(markers$gene, function(x) strsplit(x, split = "-")[[1]][1])
markers <- markers %>% left_join(xenium_lung[,c("Gene", "Annotation")], by = c("gene_symbol" = "Gene"))
topmarkers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.3, p_val_adj < 0.05, pct.1 > 0.5)

####
## annotation ####
####

# cluster 1.0
# 1 Epithelial Cells
# 2 B Cells 1
# 3 B Cells 2
# 4 Stromal Cells (C1R/C1S+)
# 5 Epithelial Cells
# 6 Stromal Cells (C1R/C1S+)
# 7 Macrophage/DC 1
# 8 T Cells 1
# 9 Endothelial Cells
# 10 Smooth muscle cell
# 11 Plasmacytoid dendritic cell
# 12 T Cells 2
# 13 unknown
# 14 Germinal Center 
# 15 Proliferating Cells 2
# 16 Macrophage/DC 2
# 17 Macrophage/DC 3
# 18 unknown
# 19 Plasma cells
# 20 T Cells 3
# 21 B Cells 3
# 22 T Cells 4
# 23 unknown
# 24 Endothelial Cells
# 25 T Cells 5 
# 26 Macrophage/DC 4
# 27 Mast Cells
# 28 Lymphatic endothelial cells
# 29 Macrophage/DC 5

clusters <- as.numeric(xenium_reg$clusters_1)
celltypes <- c(
  "Epithelial Cells",
  "B Cells 1",
  "B Cells 2",
  "Stromal Cells (C1R/C1S+)",
  "Epithelial Cells",
  "Stromal Cells (C1R/C1S+)",
  "Macrophage/DC 1",
  "T Cells 1",
  "Endothelial Cells",
  "Smooth Muscle Cells",
  "Plasmacytoid Dendritic Cells",
  "T Cells 2",
  "Unknown",
  "Germinal Center",
  "Proliferating Cells",
  "Macrophage/DC 2",
  "Macrophage/DC 3",
  "Unknown",
  "Plasma Cells",
  "T Cells 3",
  "B Cells 3",
  "T Cells 4",
  "Unknown",
  "Endothelial Cells",
  "T Cells 5",
  "Macrophage/DC 4",
  "Mast Cells",
  "Lymphatic Endothelial Cells",
  "Macrophage/DC 5"
)
xenium_reg$CellType <- celltypes[clusters]

# update germinal center
tmp <- xenium_reg$CellType
tmp[xenium_reg$Clusters_withsub == "14_1"] <- "GC Proliferating"
tmp[xenium_reg$Clusters_withsub == "14_2"] <- "GC CD21+"
xenium_reg$CellType <- tmp

# make major cell types
tmp <- xenium_reg$CellType
tmp[grepl("B Cells", tmp)] <- "B Cells"
tmp[grepl("T Cells", tmp)] <- "T Cells"
tmp[grepl("Macrophage\\/DC", tmp)] <- "Macrophage/DC"
xenium_reg$MajorCellType <- tmp

# visualize cell type
spatialpoints <- vrSpatialPoints(xenium_reg)[xenium_reg$Clusters %in% unique(xenium_reg$Clusters)[unique(xenium_reg$Clusters) != "24"]]
xenium_reg_vis <- subset(xenium_reg, spatialpoints = spatialpoints)
vrEmbeddingPlot(xenium_reg_vis, 
                group.by = "CellType",
                embedding = "umap", 
                pt.size = 0.4, label = TRUE)

# cell type spatial
colors <- hue_pal(length(unique(xenium_reg_vis$MajorCellType)))
set.seed(6)
colors <- sample(colors)
colors <- setNames(colors, unique(xenium_reg_vis$MajorCellType))
vrSpatialPlot(xenium_reg_vis, group.by = "MajorCellType", plot.segments = TRUE, alpha = 1, colors = colors)
ggsave("clustering.pdf", plot = last_plot(), device = "pdf", height = 8, width = 8)
vrSpatialPlot(xenium_reg_vis, group.by = "MajorCellType", plot.segments = TRUE, alpha = 1, colors = colors) + 
  theme(legend.position = "none") + 
  labs(title = "")
# ggsave("clustering_notitle.pdf", plot = last_plot(), device = "pdf", height = 8, width = 8)

####
### marker Analysis ####
####

# marker analysis
datax <- vrData(xenium_reg, feat_type = vrFeatureTypeNames(xenium_reg), norm = TRUE)
xenium_reg_seu <- CreateSeuratObject(as(datax, "dgCMatrix"), meta.data = as.data.frame(Metadata(xenium_reg)))
xenium_reg_seu$clusters_1 <- as.character(xenium_reg_seu$clusters_1)
new_datax <- datax
rownames(new_datax) <- rownames(xenium_reg_seu)
xenium_reg_seu <- SetAssayData(
  object = xenium_reg_seu,
  layer = "data",
  new.data = as(new_datax,"dgCMatrix"),
  assay = "RNA"
)
Idents(xenium_reg_seu) <- "CellType"
markers <- FindAllMarkers(xenium_reg_seu, features = Features(xenium_reg_seu))
markers$gene_symbol <- sapply(markers$gene, function(x) strsplit(x, split = "-")[[1]][1])
topmarkers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.3, p_val_adj < 0.05, pct.1 > 0.5)
# write.xlsx(as.data.frame(topmarkers), 
#            file = "../../Nature Methods Revision/supp table/Supplementary Table 4.xlsx", 
#            sheetName = "with_mIF")

####
## save voltron ####
####

xenium_reg <- saveVoltRon(xenium_reg)

####
## as.anndata ####
####

vrMainFeatureType(xenium_reg) <- "RNA"
xenium_reg$clusters_1 <- as.character(xenium_reg$clusters_1)
as.AnnData(xenium_reg, file = "data/Tonsil/xenium_tonsil_IF.h5ad", assay = "Xenium", flip_coordinates = TRUE)

####
## Germinal Center (Zoom-in) ####
####

# xenium_reg_sub <- subset(xenium_reg, interactive = TRUE)
# xenium_reg_sub$subset_info_list[[1]]
xenium_reg_sub <- subset(xenium_reg, image = "1687x1448+5117+9459")
xenium_reg_sub <- combineChannels(xenium_reg_sub, channels = "DAPI", colors = "white", channel_key = "DAPI_new")
xenium_reg_sub <- combineChannels(xenium_reg_sub, channels = "CD20", colors = "white", channel_key = "CD20_new")
xenium_reg_sub <- combineChannels(xenium_reg_sub, channels = "CD21", colors = "white", channel_key = "CD21_new")
xenium_reg_sub <- modulateImage(xenium_reg_sub, channel = "DAPI_new", brightness = 70)
xenium_reg_sub <- modulateImage(xenium_reg_sub, channel = "CD20_new", brightness = 700)
xenium_reg_sub <- modulateImage(xenium_reg_sub, channel = "CD21_new", brightness = 200)
xenium_reg_sub <- combineChannels(xenium_reg_sub, 
                                  channels = c("DAPI_new", "CD20_new", "CD21_new"),
                                  colors = c("blue", "red", "yellow"))
image_ggplot(vrImages(xenium_reg_sub, channel = "combined"))
vrSpatialPlot(xenium_reg_sub, group.by = "MajorCellType", plot.segments = TRUE, alpha = 0, channel = "combined") +
  labs(title = "") + 
  theme(legend.position = "none")
ggsave("multiplex.pdf", plot = last_plot(), device = "pdf", height = 4, width = 4)

# colors defined above
vrSpatialPlot(xenium_reg_sub, group.by = "MajorCellType", plot.segments = TRUE, alpha = 1,
                    group.ids = c("B Cells", "T Cells", "GC CD21+", "GC Proliferating", "Proliferating Cells"), 
              channel = "DAPI_new", background.color = "black",
              colors = colors[c("B Cells", "T Cells", "GC CD21+", "GC Proliferating", "Proliferating Cells")])+
  labs(title = "") + 
  theme(legend.position = "none")
ggsave("celltype.pdf", plot = last_plot(), device = "pdf", height = 4, width = 4)

vrMainFeatureType(xenium_reg_sub) <- "Protein"
vrSpatialFeaturePlot(xenium_reg_sub, features = c("CD21"), plot.segments = TRUE, alpha = 1, background.color = TRUE) + 
  labs(title = "") + 
  theme(legend.position = "none")
ggsave("cd21.pdf", plot = last_plot(), device = "pdf", height = 4, width = 4)

vrMainFeatureType(xenium_reg_sub) <- "RNA"
vrSpatialFeaturePlot(xenium_reg_sub, features = c("MKI67"), plot.segments = TRUE, alpha = 1, background.color = TRUE) +
  labs(title = "") + 
  theme(legend.position = "none")
ggsave("MKI67.pdf", plot = last_plot(), device = "pdf", height = 4, width = 4)

# subset 2
# xenium_reg_sub <- subset(xenium_reg, interactive = TRUE)
# xenium_reg_sub$subset_info_list[[1]]
xenium_reg_sub <- subset(xenium_reg, image = "2946x1809+4968+9163")

vrMainFeatureType(xenium_reg_sub) <- "Protein"
g1 <- vrSpatialFeaturePlot(xenium_reg_sub, features = c("CD21"), plot.segments = TRUE, alpha = 1)
vrMainFeatureType(xenium_reg_sub) <- "RNA"
g2 <- vrSpatialFeaturePlot(xenium_reg_sub, features = c("CD38", "UBE2C", "MKI67", "CXCL13"), 
                           plot.segments = TRUE, collapse.plots = FALSE, alpha = 1)
g2 <- c(g2, list(g1))
ggpubr::ggarrange(plotlist = g2, ncol = 3, nrow = 2)

# subset 3
# xenium_reg_sub <- subset(xenium_reg, interactive = TRUE)
# xenium_reg_sub$subset_info_list[[1]]
xenium_reg_sub <- subset(xenium_reg, image = "1976x2349+8532+4451")

vrMainFeatureType(xenium_reg_sub) <- "Protein"
g1 <- vrSpatialFeaturePlot(xenium_reg_sub, features = c("CD21", "CD20"), 
                           plot.segments = TRUE, alpha = 1, collapse.plots = FALSE)
vrMainFeatureType(xenium_reg_sub) <- "RNA"
g2 <- vrSpatialFeaturePlot(xenium_reg_sub, features = c("CD38", "UBE2C", "MKI67", "CXCL13"), 
                           plot.segments = TRUE, collapse.plots = FALSE, alpha = 1)
g2 <- c(g2, g1)
ggpubr::ggarrange(plotlist = g2, ncol = 3, nrow = 2)

####
# Xenium Tonsil (Adjacent Section) ####
####

####
## Import QuPath Data ####
####

library(sf)
segments <- read_sf("../data/Tonsil/IF_QuPath/measurements_featurecollection.geojson", 
                    type = 3)
segments <- generateSegments(segments)
adj_vr <- importQuPathIF(measurements = "../data/Tonsil/IF_QuPath/measurements_adjusted.txt",
                         segments = segments, 
                         image = "../../../../data/Erik/Erik_Alignment_TMA/out/Adjacent section/Core_13.ome.tif", 
                         channels = c("DAPI", "CD20", "CD68", "CD21", "FOXP3", "CD45RB"),
                         series = 1,
                         resolution = 1, 
                         sample_name = "IFdata")

####
## adjacent section alignment ####
####

# Subset Xenium data
# xen <- subset(xen, interactive = TRUE)
# xen_subset <- xen$subsets[[1]]
xen_subset <- subset(xen, image = "5149x4945+7145+4380")

# register 617.88x593.4 vs 606.95x606.8
xen_subset2 <- modulateImage(xen_subset, brightness = 500)
adj_vr2 <- modulateImage(adj_vr, brightness = 700, channel = "DAPI")
xen_reg <- registerSpatialData(object_list = list(xen_subset2, adj_vr2),
                               mapping_parameters = readRDS("../data/Tonsil/adjsection_registration_parameters.rds"),
                               interactive = TRUE)
xenium_reg <- xen_reg$registered_spat

# merge data
xenium_reg <- merge(xenium_reg[[1]], xenium_reg[[2]], samples = "XeniumBlock")

# saveRDS(xen_reg$accuracy_metrics, file = "results/adj_section_accuracy.rds")

# visualize accuracy
accuracy_metrics <- readRDS(file = "results/adj_section_accuracy.rds")
VoltRon:::.plot_ssim_map(accuracy_metrics$ssim_maps$`2`$Coarse)
VoltRon:::.plot_ssim_map(accuracy_metrics$ssim_maps$`2`$Fine)

####
## Analyze IF data ####
####

# get IF feature names
vrFeatures(xenium_reg)

# normalize Protein markers
vrMainAssay(xenium_reg) <- "IF"
xenium_reg <- normalizeData(xenium_reg, method = "hyper.arcsine", scale = 0.2)

# visualize normalized protein markers
vrSpatialFeaturePlot(xenium_reg, features = vrFeatures(xenium_reg), ncol = 3, alpha = 1)

####
## Clustering (IF) ####
####

# get umap from normalized counts
xenium_reg <- getUMAP(xenium_reg, data.type = "norm")

# visualize markers on embedding
vrEmbeddingFeaturePlot(xenium_reg, features = vrFeatures(xenium_reg), 
                       embedding = "umap" , ncol = 3)

# clustering k means
xenium_reg <- getClusters(xenium_reg, method = "kmeans", nclus = 7, label = "Clusters")
g1 <- vrEmbeddingPlot(xenium_reg, group.by = "Clusters", assay = "IF",
                embedding = "umap")
g2 <- vrSpatialPlot(xenium_reg, group.by = "Clusters", assay = "IF")
g3 <- vrEmbeddingPlot(xenium_reg, group.by = "Clusters", assay = "IF",
                      embedding = "umap", group.ids = c(4,5))
g4 <- vrSpatialPlot(xenium_reg, group.by = "Clusters", assay = "IF", group.ids = c(4,5), plot.segments = TRUE)
(g1 | g2) / (g3 | g4)

# subclustering
xenium_reg_sub <- subset(xenium_reg, subset = Clusters %in% c(1,2,3,5,6))
xenium_reg_sub <- getClusters(xenium_reg_sub, method = "kmeans", nclus = 7, label = "Clusters")
vrEmbeddingPlot(xenium_reg_sub, group.by = "Clusters", 
                embedding = "umap")
xenium_reg_sub <- subset(xenium_reg_sub, subset = Clusters %in% c(3,5,6,7))
vrEmbeddingPlot(xenium_reg_sub, group.by = "Clusters", 
                embedding = "umap")
xenium_reg_sub <- getClusters(xenium_reg_sub, method = "kmeans", nclus = 8, label = "subclusters")
vrEmbeddingPlot(xenium_reg_sub, group.by = "subclusters",
                embedding = "umap")
xenium_reg_sub <- subset(xenium_reg_sub, subset = subclusters %in% c(1,3,4,5,6,7,8))
xenium_reg_sub <- getClusters(xenium_reg_sub, method = "kmeans", nclus = 7, label = "Clusters")
g1 <- vrEmbeddingPlot(xenium_reg_sub, group.by = "Clusters", 
                embedding = "umap")
g2 <- vrEmbeddingFeaturePlot(xenium_reg_sub, 
                       features = c("CD21", "CD20", "FOXP3", "CD45RB"), 
                       embedding = "umap" , ncol = 2)
g1 | g2
  
# insert clusters
clusters <- setNames(rep("Other", length(vrSpatialPoints(xenium_reg))),
                     vrSpatialPoints(xenium_reg))
clusters[xenium_reg$Clusters == 4] <- "B Cells"
clusters[vrSpatialPoints(xenium_reg_sub)] <- xenium_reg_sub$Clusters
xenium_reg$annotation <- clusters

# cell type annotation
vrMainAssay(xenium_reg) <- "IF"
celltype <- xenium_reg$annotation 
celltype[celltype == 5] <- "Follicular Dendritic Cells"
celltype[celltype == 7] <- "FOXP3+"
celltype[celltype %in% c(2,3,6)] <- "B Cells"
celltype[celltype == 4] <- "GC (CD45RB-)"
celltype[celltype == 1] <- "GC (CD45RB+)"
xenium_reg$CellType <- celltype
celltype <- xenium_reg_sub$Clusters 
celltype[celltype == 5] <- "Follicular Dendritic Cells"
celltype[celltype == 7] <- "FOXP3+"
celltype[celltype %in% c(2,3,6)] <- "B Cells"
celltype[celltype == 4] <- "GC (CD45RB-)"
celltype[celltype == 1] <- "GC (CD45RB+)"
xenium_reg_sub$CellType <- celltype

# visualize
colors <- hue_pal(length(unique(xenium_reg$CellType)))
names(colors) <- c(unique(xenium_reg$CellType))
vrEmbeddingPlot(xenium_reg_sub, group.by = "CellType", 
                      embedding = "umap", colors = colors[names(colors) %in% unique(xenium_reg_sub$CellType)])
ggsave("protein_embedding.pdf", plot = last_plot(), device = "pdf", height = 8, width = 10)
vrSpatialPlot(xenium_reg, group.by = "CellType",
              plot.segments = FALSE, alpha = 1)
ggsave("protein_clustering.pdf", plot = last_plot(), device = "pdf", height = 8, width = 8)
vrSpatialPlot(xenium_reg, group.by = "CellType",
              plot.segments = FALSE, alpha = 1) + 
  theme(legend.position = "none") + 
  labs(title = "")
ggsave("protein_clustering_notitle.pdf", plot = last_plot(), device = "pdf", height = 8, width = 8)

# heatmap visualization
# vrHeatmapPlot(xenium_reg, features = c("CD21", "CD20", "FOXP3", "CD45RB"), group.by = "CellType")

# saveRDS(xenium_reg, "../data/Tonsil/xenium_reg.rds")
xenium_reg <- readRDS("../data/Tonsil/xenium_reg.rds")

# visualize non-others 
xenium_reg_sub <- subset(xenium_reg, subset = CellType %in% unique(xenium_reg$CellType)[unique(xenium_reg$CellType) != "Other"])
g1 <- vrEmbeddingPlot(xenium_reg_sub, group.by = "CellType", 
                      embedding = "umap")
g2 <- vrEmbeddingFeaturePlot(xenium_reg_sub, 
                             features = c("CD21", "CD20", "FOXP3", "CD45RB"), 
                             embedding = "umap" , ncol = 2)
g3 <- vrSpatialPlot(xenium_reg, group.by = "CellType", n.tile = 300)
g1 | g2

####
## Load data ####
####

xenium_reg <- readRDS("../data/Tonsil/xenium_reg.rds")
xen_if <- loadVoltRon("../data/Tonsil/Xenium_IF/")

####
## Neighborhood Analysis ####
####

# get Xenium clusters
xen_celltype <- setNames(xen_if$MajorCellType, vrSpatialPoints(xen_if))
xen_celltype[grepl("^GC",xen_celltype)] <- "Proliferating Cells (CD38+)"

# create new column
vrMainAssay(xenium_reg) <- "Xenium"
celltypes <- setNames(xenium_reg$CellType,
                      vrSpatialPoints(xenium_reg))
celltypes[names(xen_celltype)] <- as.character(xen_celltype)
celltypes[celltypes == ""] <- "Unknown"
xenium_reg$CellType <- celltypes

# visualize cell types of both groups
g1 <- vrSpatialPlot(xenium_reg, assay = "IF", group.by = "CellType")
g2 <- vrSpatialPlot(xenium_reg, assay = "Xenium", group.by = "CellType")c("GC (CD45RB-)", "GC (CD45RB+)")

# visualize Germinal Center of both groups
g1 <- vrSpatialPlot(xenium_reg, assay = "IF", group.by = "CellType", group.ids = c("GC (CD45RB-)", "GC (CD45RB+)"), n.tile = 200)
g2 <- vrSpatialPlot(xenium_reg, assay = "Xenium", group.by = "CellType", group.ids = c("Proliferating Cells (CD38+)"), n.tile = 200)
g1 | g2

# attach connectivity 
# micron to pixel ratio is 0.85
xenium_reg2 <- xenium_reg
xenium_reg2[["XeniumBlock"]]@zlocation[2] <- 5/0.85
xenium_reg2[["XeniumBlock"]]@adjacency[1,2] <- 
  xenium_reg2[["XeniumBlock"]]@adjacency[2,1] <- 1

# spatial proximity
xenium_reg2 <- getSpatialNeighbors(xenium_reg2, assay = c("Assay1", "Assay2"), radius = 20, method = "radius")
vrGraphNames(xenium_reg2)

# neighborhood analysis
results <- VoltRon::vrNeighbourhoodEnrichment(xenium_reg2, group.by = list(Assay1 = "CellType", 
                                                                           Assay2 = "CellType"), assay = c("Assay1", "Assay2"), 
                                              num.sim = 100, graph.type = "radius")
results_new <- results[results$from_value != "Other" & results$to_value != "Other",]
results_new <- results_new[results_new$from_value != "Unknown" & results_new$to_value != "Unknown",]
pdf("neighbourhood_enrichment.pdf", width = 12, height = 10)
VoltRon::vrNeighbourhoodEnrichmentPlot(results_new, assay = c("Assay1,Assay2"))
dev.off()

####
## Germinal Center (Zoom-in) ####
####

# colors
spatialpoints <- vrSpatialPoints(xenium_reg)[xenium_reg$Clusters %in% unique(xenium_reg$Clusters)[unique(xenium_reg$Clusters) != "24"]]
xenium_reg_vis <- subset(xenium_reg, spatialpoints = spatialpoints)
colors <- hue_pal(length(unique(xenium_reg_vis$CellType)))
set.seed(6)
colors <- sample(colors)
colors <- setNames(colors, unique(xenium_reg_vis$CellType))

xenium_reg_sub <- subset(xenium_reg, assay = "Xenium")
# xenium_reg_sub <- subset(xenium_reg_sub, interactive = TRUE)
# xenium_reg_sub$subset_info_list[[1]]
xenium_reg_sub <- subset(xenium_reg_sub, image = "712x707+3713+2087")
colors_vis <- colors[unique(xenium_reg_sub$CellType)]
g1 <- vrSpatialPlot(xenium_reg_sub, plot.segments = TRUE, group.by = "CellType", colors = colors_vis, background.color = "black") 
ggsave("xenium_zoom.pdf", plot = g1, device = "pdf", height = 4, width = 6)
g1 <- g1 + labs(title = "") + 
  theme(legend.position = "none")
ggsave("xenium_zoom_nolabel.pdf", plot = g1, device = "pdf", height = 4, width = 4)

xenium_reg_sub <- subset(xenium_reg, assay = "IF")
xenium_reg_sub <- subset(xenium_reg_sub, image = "712x707+3713+2087")
g1 <- vrSpatialPlot(xenium_reg_sub, plot.segments = TRUE, group.by = "CellType", background.color = "black")
ggsave("if_zoom.pdf", plot = g1, device = "pdf", height = 4, width = 6)
g1 <- g1 + labs(title = "") + 
  theme(legend.position = "none")
ggsave("if_zoom_nolabel.pdf", plot = g1, device = "pdf", height = 4, width = 4)

####
## marker analysis ####
####

# marker analysis
datax <- vrData(xenium_reg, feat_type = vrFeatureTypeNames(xenium_reg), norm = TRUE)
xenium_reg_seu <- CreateSeuratObject(as(datax, "dgCMatrix"), meta.data = as.data.frame(Metadata(xenium_reg)))
new_datax <- datax
rownames(new_datax) <- rownames(xenium_reg_seu)
xenium_reg_seu <- SetAssayData(
  object = xenium_reg_seu,
  layer = "data",
  new.data = as(new_datax,"dgCMatrix"),
  assay = "RNA"
)
Idents(xenium_reg_seu) <- "CellType"
markers <- FindAllMarkers(xenium_reg_seu, features = Features(xenium_reg_seu))
markers$gene_symbol <- sapply(markers$gene, function(x) strsplit(x, split = "-")[[1]][1])
topmarkers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.3, p_val_adj < 0.05, pct.1 > 0.5)
# write.xlsx(as.data.frame(topmarkers), 
#            file = "../../Nature Methods Revision/supp table/Supplementary Table 4.xlsx", 
#            sheetName = "without_mIF", append = TRUE)