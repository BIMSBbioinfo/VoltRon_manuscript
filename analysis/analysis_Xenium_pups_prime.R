####
# Xenium 5000 ####
####

####
## Import Xenium ####
####


Xen_R1 <- importXenium("../../../../data/xenium/Xenium_V1_mouse_pup_outs/", sample_name = "XeniumR1",
                       resolution_level = 3, overwrite_resolution = TRUE, import_molecules = FALSE)

####
## save to disk ####
####

# load other libraries
library(BPCells)
library(ImageArray)

# save to disk
Xen_R1_disk <- saveVoltRon(Xen_R1, format = "HDF5VoltRon", output = "../data/Xenium_pups_prime", replace = TRUE)
Xen_R1_disk <- loadVoltRon("../data/Xenium_pups_prime//")

####
## Processing ####
####

# filter
spatialpoints <- vrSpatialPoints(Xen_R1_disk)[as.vector(Metadata(Xen_R1_disk)$Count > 100)]
Xen_R1_disk <- subset(Xen_R1_disk, spatialpoints = spatialpoints)

# normalize
Xen_R1_disk <- normalizeData(Xen_R1_disk, sizefactor = 10000)

# Visualize skin, heart and brain
vrSpatialFeaturePlot(Xen_R1_disk, features = c("Dsc3", "Gfap"), n.tile = 400, alpha = 1, ncol = 3, log = TRUE)

# select features
Xen_R1_disk <- getFeatures(Xen_R1_disk, n = 3000)
selectedfeatures <- getVariableFeatures(Xen_R1_disk)

# PCA and UMAP
Xen_R1_disk <- getPCA(Xen_R1_disk, dims = 30, overwrite = TRUE, features = selectedfeatures)
Xen_R1_disk <- getUMAP(Xen_R1_disk, dims = 1:30)
# saveVoltRon(Xen_R1_disk)

####
## clustering ####
####

# clustering
Xen_R1_disk <- getProfileNeighbors(Xen_R1_disk, dims = 1:30, k = 10, method = "SNN")
Xen_R1_disk <- getClusters(Xen_R1_disk, resolution = 0.5, label = "Clusters", graph = "SNN")

# visualize
library(patchwork)
g1 <- vrSpatialPlot(Xen_R1_disk, group.by = "Clusters", n.tile = 400, alpha = 1, legend.loc = "none")
g2 <- vrEmbeddingPlot(Xen_R1_disk, group.by = "Clusters", embedding = "umap", label = TRUE)

####
## marker analysis ####
####

# marker analysis
xenium_markers <- read.csv("../../../test/mainWorkflows/data/Xenium_mMulti_v1_metadata_annotations.csv")
xenium_markers2 <- read.csv("../../../test/mainWorkflows/data/XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv")
Xen_R1_disk_seu <- as.Seurat(Xen_R1_disk, cell.assay = "Xenium", type = "image")
Idents(Xen_R1_disk_seu) <- "Clusters"
Xen_R1_disk_seu <- NormalizeData(Xen_R1_disk_seu, scale.factor = 10000)
markers <- Seurat::FindAllMarkers(Xen_R1_disk_seu)
# saveRDS(markers, file = "data/xeniumpupsprime_markers.rds")
markers <- readRDS("data/xeniumpupsprime_markers.rds")
markers <- markers %>% left_join(xenium_markers[,c("Gene", "Tissues", "Cell.types")], by = c("gene" = "Gene"))
markers <- markers %>% left_join(xenium_markers2[,c("gene_name", "cell_type")], by = c("gene" = "gene_name"))
topmarkers <- markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1, p_val_adj < 0.05, pct.1 > 0.5) %>% 
  top_n(n = 20, wt = avg_log2FC)

####
## annotation ####
####

Xen_R1_disk2 <- loadVoltRon("data/Xenium_pups_prime/")
vrSpatialPlot(Xen_R1_disk, group.by = "Clusters", n.tile = 400, alpha = 1, legend.loc = "none", group.ids = "1")
# cluster 1 - Skin (Epidermis)
# cluster 2 - Adipose
# cluster 3 - Immune 1
# cluster 4 - Skin (Dermis)
# cluster 5 - Other
# cluster 6 - Other
# cluster 7 - Epithelium
# cluster 8 - Muscle (Cardiac/Skeletal)
# cluster 9 - Central Nervous System
# cluster 10 - Bone/connective tissue ECM  
# cluster 11 - Brain
# cluster 12 - Neural progenitor cells (NPCs)
# cluster 13 - Vasculature (Smooth muscle)
# cluster 14 - Kidney/Intestine
# cluster 15 - Liver
# cluster 16 - Cartilage
# cluster 17 - Vasculature (Lymphatic)
# cluster 18 - Immune 2
# cluster 19 - Blood/Hematopoietic
# cluster 20 - Vasculature (Endothelium)
# cluster 21 - Eye (Lens & Retina) 
# cluster 22 - Lung
annotation <- c(
  "Skin (Epidermis)",
  "Adipose Tissue",
  "Immune 1",
  "Skin (Dermis)",
  "Other",
  "Other",
  "Epithelium",
  "Muscle (Cardiac/Skeletal)",
  "Central Nervous System",
  "Bone/connective tissue ECM",
  "Brain",
  "Neural progenitor cells (NPCs)",
  "Vasculature (Smooth muscle)",
  "Kidney/Intestine",
  "Liver",
  "Cartilage",
  "Vasculature (Lymphatic)",
  "Immune 2",
  "Blood/Hematopoietic",
  "Vasculature (Endothelium)",
  "Eye (Lens & Retina)",
  "Lung"
)
Xen_R1_disk$CellType <- annotation[as.numeric(Xen_R1_disk$Clusters)]

vrSpatialPlot(Xen_R1_disk, group.by = "CellType", n.tile = 400, alpha = 1, legend.loc = "none")