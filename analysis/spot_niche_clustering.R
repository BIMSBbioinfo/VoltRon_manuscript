# library(VoltRon)
library(philentropy)
library(patchwork)
library(ggpubr)
library(mclust)
library(Seurat)
library(dplyr)

####
## Import data ####
####

DLPFC_merged <- readRDS("../data/DLFPC/Visium&Visium_data_decon_registered.rds")

####
## get ground truth ####
####

DLPFC_1_labels <- readRDS("../data/DLFPC/151673_labels.rds")
DLPFC_2_labels <- readRDS("../data/DLFPC/151674_labels.rds")
DLPFC_3_labels <- readRDS("../data/DLFPC/151675_labels.rds")
DLPFC_4_labels <- readRDS("../data/DLFPC/151676_labels.rds")
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

# Visualize
vrMainFeatureType(DLPFC_merged) <- "Decon"
vrSpatialFeaturePlot(DLPFC_merged, assay = c("Assay1", "Assay4"), 
                     features = c("Excit_E", "Excit_A", "Inhib_C"),
                     crop = TRUE, ncol = 3, alpha = 1, keep.scale = "all")

####
## 2D Clustering ####
####

# DLPFC_merged 2d
DLPFC_merged_2d <- DLPFC_merged

# remove adjacency in the tissue block
tmp <- DLPFC_merged_2d@samples$DLPFC_Block@adjacency
tmp[1:4,1:4] <- diag(4)
DLPFC_merged_2d@samples$DLPFC_Block@adjacency <- tmp

####
### Niche Assay from Decon ####
####

DLPFC_merged_2d <- getSpatialNeighbors(DLPFC_merged_2d, method = "radius")
vrMainFeatureType(DLPFC_merged_2d) <- "Decon"
DLPFC_merged_2d <- getNicheAssay(DLPFC_merged_2d, graph.type = "radius")
vrMainFeatureType(DLPFC_merged_2d) <- "Niche"

####
### Processing ####
####

vrMainFeatureType(DLPFC_merged_2d) <- "Niche"
DLPFC_merged_2d <- normalizeData(DLPFC_merged_2d, method = "CLR")

####
### Clustering ####
####

# embedding
DLPFC_merged_2d <- getUMAP(DLPFC_merged_2d, data.type = "norm")
vrEmbeddingPlot(DLPFC_merged_2d, embedding = "umap", group.by = "Sample")

# clustering 
DLPFC_merged_2d <- getProfileNeighbors(DLPFC_merged_2d, data.type = "norm", method = "SNN")
DLPFC_merged_2d <- getClusters(DLPFC_merged_2d, resolution = 0.49, graph = "SNN", label = "clusters_SNN")

# clustering K Means
DLPFC_merged_2d <- getClusters(DLPFC_merged_2d, method = "kmeans", nclus = 7, label = "clusters_kmeans")

# clustering Manhattan
DLPFC_merged_2d <- getClusters(DLPFC_merged_2d, method = "hierarchical", nclus = 7, distance_measure = "manhattan", label = "clusters_hier")

# clustering JSD
vrdata <- t(vrData(DLPFC_merged_2d, norm = FALSE))
propor_dis <- philentropy::distance(vrdata, method = "jensen-shannon")
rownames(propor_dis) <- colnames(propor_dis) <- rownames(vrdata)
propor_dis <- as.dist(propor_dis)
clusters <- stats::hclust(d = propor_dis, method = "ward.D2")
clusters <- stats::cutree(clusters, k = 7)
clusters <- list(names = names(clusters), membership = clusters)
spatialpoints <- vrSpatialPoints(DLPFC_merged_2d)
membership <- setNames(rep(NA,length(spatialpoints)), spatialpoints)
membership[clusters$names] <- clusters$membership
DLPFC_merged_2d <- addMetadata(DLPFC_merged_2d, value = membership, label = "clusters_hierjsd")

# save clustered
# saveRDS(DLPFC_merged_2d, file = "../data/DLFPC/DLPFC_merged_nicheclustered_2d.rds")
# DLPFC_merged_2d <- readRDS("../data/DLFPC/DLPFC_merged_nicheclustered_2d.rds")

# visualize
colors <- hue_pal(7)
names(colors) <- c(5,7,6,2,4,1,3)
vrSpatialPlot(DLPFC_merged_2d, group.by = "clusters_kmeans", alpha = 1, nrow = 2, crop = TRUE, colors = colors)
ggsave(filename = "../../Nature Methods Revision/Images/integration/images/registeration/registeration_visium/DLPFC_nicheclusters_2d.pdf", 
       plot = last_plot(), device = "pdf", width = 10, height = 8, units = "in")

# 

####
## 3D Clustering ####
####

####
### Niche Assay from Decon ####
####

DLPFC_merged <- getSpatialNeighbors(DLPFC_merged, method = "radius")
vrMainFeatureType(DLPFC_merged) <- "Decon"
DLPFC_merged <- getNicheAssay(DLPFC_merged, graph.type = "radius")
vrMainFeatureType(DLPFC_merged) <- "Niche"

####
### Processing ####
####

vrMainFeatureType(DLPFC_merged) <- "Niche"
DLPFC_merged <- normalizeData(DLPFC_merged, method = "CLR")

####
### Clustering ####
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
# saveRDS(DLPFC_merged, file = "../data/DLFPC/DLPFC_merged_nicheclustered.rds")
# DLPFC_merged <- readRDS("../data/DLFPC/DLPFC_merged_nicheclustered.rds")

# colors
colors_layer <- hue_pal(8)
names(colors_layer) <- unique(DLPFC_merged$Layers)

# visualize assay 1
g_list <- list()
g_list[[1]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay1", group.by = "clusters_SNN", ncol = 1)
g_list[[2]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay1", group.by = "clusters_kmeans", ncol = 1)
g_list[[3]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay1", group.by = "clusters_hier", ncol = 1)
g_list[[4]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay1", group.by = "clusters_hierjsd", ncol = 1)
g_list[[5]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay1", group.by = "Layers", ncol = 1, colors = colors_layer)
ggpubr::ggarrange(plotlist = g_list, ncol = 3, nrow = 2)

# visualize assay 4
g_list2 <- list()
g_list2[[1]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay3", group.by = "clusters_SNN", ncol = 1)
g_list2[[2]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay3", group.by = "clusters_kmeans", ncol = 1)
g_list2[[3]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay3", group.by = "clusters_hier", ncol = 1)
g_list2[[4]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay3", group.by = "clusters_hierjsd", ncol = 1)
g_list2[[5]] <- vrSpatialPlot(DLPFC_merged, assay = "Assay3", group.by = "Layers", ncol = 1, colors = colors_layer)
ggpubr::ggarrange(plotlist = g_list2, ncol = 3, nrow = 2)

g_list_all <- c(g_list, g_list2)
ggpubr::ggarrange(plotlist = g_list_all, ncol = 5, nrow = 2)
# ggsave(filename = "../../Nature Methods Revision/Images/Supplementary Material/SpatiallyAwareAnalysis/spot_comparison.pdf", 
#        plot = last_plot(), device = "pdf", width = 30, height = 10, units = "in")

####
### Heatmap ####
####

# visualize heatmap
g1 <- vrHeatmapPlot(DLPFC_merged, features = vrFeatures(DLPFC_merged), group.by = "clusters_kmeans")
g2 <- vrHeatmapPlot(DLPFC_merged, features = vrFeatures(DLPFC_merged), group.by = "clusters_hierjsd")
g3 <- vrHeatmapPlot(DLPFC_merged, features = vrFeatures(DLPFC_merged), group.by = "Layers")
g1 + g2 + g3

####
## ARI ####
####

DLPFC_merged_2d <- readRDS("../data/DLFPC/DLPFC_merged_nicheclustered_2d.rds")
DLPFC_merged <- readRDS("../data/DLFPC/DLPFC_merged_nicheclustered.rds")

####
### 2D ####
####

sample_metadata <- SampleMetadata(DLPFC_merged_2d)
sample_metadata <- data.frame(sample_metadata, labels = c("151673", "151674", "151675", "151676"))
metadata <- Metadata(DLPFC_merged_2d)
variables <- colnames(metadata)
variables <- variables[grepl("Layers|^clusters", variables)]
metadata <- metadata[metadata$Layers != "none",]

results_list_2d <- list()
for(assy in unique(metadata$assay_id)){
  cur_metadata <- metadata[metadata$assay_id == assy,]
  results <- matrix(1, nrow = length(variables), ncol = length(variables))
  for(i in 1:(length(variables)-1)){
    for(j in (i+1):length(variables)){
      results[i,j] <- results[j,i] <- 
        mclust::adjustedRandIndex(cur_metadata[[variables[i]]], cur_metadata[[variables[j]]])
    }
  }
  rownames(results) <- colnames(results) <- variables
  results_list_2d[[sample_metadata[assy, "labels"]]] <- results[-1,"Layers"]
}
results_list_2d <- do.call(results_list_2d, what = "rbind")

####
### 3D ####
####

sample_metadata <- SampleMetadata(DLPFC_merged)
sample_metadata <- data.frame(sample_metadata, labels = c("151673", "151674", "151675", "151676"))
metadata <- Metadata(DLPFC_merged)
variables <- colnames(metadata)
variables <- variables[grepl("Layers|^clusters", variables)]
metadata <- metadata[metadata$Layers != "none",]

results_list <- list()
for(assy in unique(metadata$assay_id)){
  cur_metadata <- metadata[metadata$assay_id == assy,]
  results <- matrix(1, nrow = length(variables), ncol = length(variables))
  for(i in 1:(length(variables)-1)){
    print(variables[i])
    for(j in (i+1):length(variables)){
      print(c(i,j))
      results[i,j] <- results[j,i] <- 
        mclust::adjustedRandIndex(cur_metadata[[variables[i]]], cur_metadata[[variables[j]]])
    }
  }
  rownames(results) <- colnames(results) <- variables
  print(results)
  results_list[[sample_metadata[assy, "labels"]]] <- results[-1,"Layers"]
}
results_list <- do.call(results_list, what = "rbind")

####
### visualize ####
####

# visualize
results_list <- reshape2::melt(results_list)
colnames(results_list) <- c("Sample", "Method", "ARI")
results_list_2d <- reshape2::melt(results_list_2d)
colnames(results_list_2d) <- c("Sample", "Method", "ARI")
results_list_merged <- data.frame(rbind(results_list_2d, results_list), 
                                  Type = c(rep("2D", nrow(results_list_2d)), rep("3D", nrow(results_list))))
tmp <- results_list_merged$Sample
tmp <- factor(tmp, levels = c("151673", "151674", "151675", "151676"))
results_list_merged$Sample <- tmp
ggplot(results_list_merged, aes(x = Method, y = ARI, fill = Sample)) + 
  geom_bar(size = 5, stat = "identity", position = position_dodge()) + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust = 1)) +
  ylab("") + xlab("")+
  facet_grid(. ~ Type) +
  theme_classic() + 
  ylim(0,0.8) + 
  theme(axis.text.x = element_text(size=7, angle=45, hjust=1, vjust = 1),
        axis.text.y = element_text(size=7)) +
  scale_fill_manual(values = c("#440154", "#21908C", "#FDE725", "purple")) 
# ggsave(filename = "../../Nature Methods Revision/Images/Supplementary Material/SpatiallyAwareAnalysis/spot_comparison_ARI.pdf", 
#        plot = last_plot(), device = "pdf", width = 8, height = 4, units = "in")

# test results
# write.table(results_list_merged, file = "../data/DLFPC/DLPFC_merged_nicheclustered_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
results_list_merged <- read.table("../data/DLFPC/DLPFC_merged_nicheclustered_results.tsv", header = TRUE)

####
## test ####
####

library(glmmTMB)
library(emmeans)
m <- glmmTMB(ARI ~ Method * Type + (1 | Sample),
             family = beta_family(link = "logit"),
             data = results_list_merged)
emmeans(m, ~ Type | Method, lmer.df = c("kenward-roger")) |> contrast("revpairwise") |> summary(adjust = "holm")