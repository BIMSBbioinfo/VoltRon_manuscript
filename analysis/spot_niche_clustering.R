library(VoltRon)
library(philentropy)
library(patchwork)
library(ggpubr)
library(mclust)
library(Seurat)
library(dplyr)

####
## Import data ####
####

DLPFC_1 <- importVisium("../../../data/10X_Visium_DLPFC/151673", sample_name = "DLPFC_1")
DLPFC_2 <- importVisium("../../../data/10X_Visium_DLPFC/151674", sample_name = "DLPFC_2")
DLPFC_3 <- importVisium("../../../data/10X_Visium_DLPFC/151675", sample_name = "DLPFC_3")
DLPFC_4 <- importVisium("../../../data/10X_Visium_DLPFC/151676", sample_name = "DLPFC_4")

DLPFC_list <- list(DLPFC_1, DLPFC_2, DLPFC_3, DLPFC_4)
DLPFC_merged <- merge(DLPFC_list[[1]], DLPFC_list[-1])

####
## get ground truth ####
####

DLPFC_1_labels <- readRDS("data/DLFPC/151673_labels.rds")
DLPFC_2_labels <- readRDS("data/DLFPC/151673_labels.rds")
DLPFC_3_labels <- readRDS("data/DLFPC/151673_labels.rds")
DLPFC_4_labels <- readRDS("data/DLFPC/151673_labels.rds")
names(DLPFC_1_labels) <- paste0(names(DLPFC_1_labels) , "_Assay1")
names(DLPFC_2_labels) <- paste0(names(DLPFC_2_labels) , "_Assay2")
names(DLPFC_3_labels) <- paste0(names(DLPFC_3_labels) , "_Assay3")
names(DLPFC_4_labels) <- paste0(names(DLPFC_4_labels) , "_Assay4")
DLPFC_labels <- c(DLPFC_1_labels,
                  DLPFC_2_labels,
                  DLPFC_3_labels,
                  DLPFC_4_labels)
tmp <- as.character(DLPFC_labels[vrSpatialPoints(DLPFC_merged)])
tmp[is.na(tmp)] <- "none"
DLPFC_merged$Layers <- tmp

# visualize
vrSpatialPlot(DLPFC_merged, group.by = "Layers", ncol = 2)  

####
## get reference ####
####

load("../../../data/10X_Visium_DLPFC/SCE_DLPFC-n3_tran-etal.rda")

####
## Deconvolution ####
####

# prepare reference
tab <- table(sce.dlpfc.tran$cellType)
sce.dlpfc.tran <- sce.dlpfc.tran[,sce.dlpfc.tran$cellType %in% names(tab)[tab > 25]]

# Deconvolute Visium spots
library(spacexr)
DLPFC_merged <- getDeconvolution(DLPFC_merged, sc.object = sce.dlpfc.tran, sc.cluster = "cellType", max_cores = 6)
# saveRDS(DLPFC_merged, file = "data/DLFPC/DLPFC_merged_decon.rds")
DLPFC_merged <- readRDS("data/DLFPC/DLPFC_merged_decon.rds")

# Visualize
vrMainFeatureType(DLPFC_merged) <- "Decon"
vrSpatialFeaturePlot(DLPFC_merged, features = c("Astro", "OPC", "Excit_E", "Excit_A"),
                     crop = TRUE, ncol = 4, alpha = 1, keep.scale = "all")

####
## Niche Assay from Decon ####
####

DLPFC_merged <- getSpatialNeighbors(DLPFC_merged, method = "radius")
vrMainFeatureType(DLPFC_merged) <- "Decon"
DLPFC_merged <- getNicheAssay(DLPFC_merged, graph.type = "radius")

####
## Processing ####
####

vrMainFeatureType(DLPFC_merged) <- "Niche"
DLPFC_merged <- normalizeData(DLPFC_merged, method = "CLR")

####
## Visualization ####
####

# embedding
DLPFC_merged <- getUMAP(DLPFC_merged, data.type = "norm")
vrEmbeddingPlot(DLPFC_merged, embedding = "umap", group.by = "Sample")

# clustering 
DLPFC_merged <- getProfileNeighbors(DLPFC_merged, data.type = "norm", method = "SNN")
DLPFC_merged <- getClusters(DLPFC_merged, resolution = 0.49, graph = "SNN", label = "clusters_SNN")

# clustering K Means
DLPFC_merged <- getClusters(DLPFC_merged, method = "kmeans", nclus = 7, label = "clusters_kmeans")

# clustering Manhattan
DLPFC_merged <- getClusters(DLPFC_merged, method = "hierarchical", nclus = 7, distance_measure = "manhattan", label = "clusters_hier")

# clustering JSD
vrdata <- t(vrData(DLPFC_merged, norm = FALSE))
propor_dis <- philentropy::distance(vrdata, method = "jensen-shannon")
rownames(propor_dis) <- colnames(propor_dis) <- rownames(vrdata)
propor_dis <- as.dist(propor_dis)
clusters <- stats::hclust(d = propor_dis, method = "ward.D2")
clusters <- stats::cutree(clusters, k = 7)
clusters <- list(names = names(clusters), membership = clusters)
spatialpoints <- vrSpatialPoints(DLPFC_merged)
membership <- setNames(rep(NA,length(spatialpoints)), spatialpoints)
membership[clusters$names] <- clusters$membership
DLPFC_merged <- addMetadata(DLPFC_merged, value = membership, label = "clusters_hierjsd")

# save clustered
saveRDS(DLPFC_merged, file = "data/DLFPC/DLPFC_merged_nicheclustered.rds")

# visualize assay 1
g_list <- list()
g_list[[1]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay1", group.by = "clusters_SNN", ncol = 1)
g_list[[2]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay1", group.by = "clusters_kmeans", ncol = 1)
g_list[[3]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay1", group.by = "clusters_hier", ncol = 1)
g_list[[4]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay1", group.by = "clusters_hierjsd", ncol = 1)
g_list[[5]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay1", group.by = "Layers", ncol = 1)
ggpubr::ggarrange(plotlist = g_list, ncol = 3, nrow = 2)

# visualize assay 1
g_list <- list()
g_list[[1]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay4", group.by = "clusters_SNN", ncol = 1)
g_list[[2]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay4", group.by = "clusters_kmeans", ncol = 1)
g_list[[3]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay4", group.by = "clusters_hier", ncol = 1)
g_list[[4]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay4", group.by = "clusters_hierjsd", ncol = 1)
g_list[[5]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay4", group.by = "Layers", ncol = 1)
ggpubr::ggarrange(plotlist = g_list, ncol = 3, nrow = 2)

####
## Heatmap ####
####

# visualize heatmap
g1 <- vrHeatmapPlot(DLPFC_merged, features = vrFeatures(DLPFC_merged), group.by = "clusters_kmeans")
g2 <- vrHeatmapPlot(DLPFC_merged, features = vrFeatures(DLPFC_merged), group.by = "clusters_hierjsd")
g3 <- vrHeatmapPlot(DLPFC_merged, features = vrFeatures(DLPFC_merged), group.by = "Layers")
g1 + g2 + g3

####
## ARI ####
####

metadata <- Metadata(DLPFC_merged)
variables <- colnames(metadata)
variables <- variables[grepl("Layers|^clusters", variables)]
metadata <- metadata[metadata$Layers != "none",]

results_list <- list()
for(assy in unique(metadata$assay_id)){
  cur_metadata <- metadata[metadata$assay_id == assy,]
  results <- matrix(1, nrow = length(variables), ncol = length(variables))
  for(i in 1:(length(variables)-1)){
    for(j in (i+1):length(variables)){
      results[i,j] <- results[j,i] <- mclust::adjustedRandIndex(cur_metadata[[variables[i]]], cur_metadata[[variables[j]]])
    }
  }
  rownames(results) <- colnames(results) <- variables
  results_list[[assy]] <- results
}
