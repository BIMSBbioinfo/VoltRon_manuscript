library(VoltRon)
library(reshape2)
library(ggplot2)

###
# Xenium 5000 panel ####
###

# import
Xen_Pros <- importXenium("../../../../data/xenium/Xenium_Prime_Human_Prostate_FFPE_outs/", sample_name = "XeniumR1",
                       resolution_level = 3, overwrite_resolution = TRUE, import_molecules = FALSE)

###
## save to disk ####
###

# time
savedisk.time <- list()

# save with BPCells
start <- proc.time()
Xen_Pros_bpcells <- saveVoltRon(Xen_Pros,
                              output = "data/Xen_Pros_bpcells/",
                              format = "HDF5VoltRon",
                              replace = TRUE,
                              verbose = FALSE)
savedisk.time[["BPCells"]] <- (proc.time() - start)[3]

# save with HDF5 DelayedArray
start <- proc.time()
Xen_Pros_hdf5 <- saveVoltRon(Xen_Pros,
                           output = "data/Xen_Pros_hdf5/",
                           format = "HDF5VoltRon",
                           replace = TRUE,
                           verbose = FALSE, 
                           as.sparse = TRUE,
                           feature.vs.obs.engine = "DelayedArray")
savedisk.time[["HDF5"]] <- (proc.time() - start)[3]

# save with HDF5 DelayedArray
start <- proc.time()
Xen_Pros_hdf5_dense <- saveVoltRon(Xen_Pros,
                                 output = "data/Xen_Pros_hdf5_dense/",
                                 format = "HDF5VoltRon",
                                 replace = TRUE,
                                 verbose = FALSE, 
                                 as.sparse = FALSE,
                                 feature.vs.obs.engine = "DelayedArray")
savedisk.time[["HDF5_dense"]] <- (proc.time() - start)[3]

# save with Zarr
start <- proc.time()
Xen_Pros_zarr <- saveVoltRon(Xen_Pros,
                           output = "data/Xen_Pros_zarr/",
                           format = "ZarrVoltRon",
                           replace = TRUE,
                           verbose = FALSE)
savedisk.time[["Zarr"]] <- (proc.time() - start)[3]

###
## workflow ####
###

# load
Xen_Pros_bpcells <- loadVoltRon(dir = "data/Xen_Pros_bpcells/")
Xen_Pros_hdf5 <- loadVoltRon(dir = "data/Xen_Pros_hdf5/")
Xen_Pros_hdf5_dense <- loadVoltRon(dir = "data/Xen_Pros_hdf5_dense/")
Xen_Pros_zarr <- loadVoltRon(dir = "data/Xen_Pros_zarr/")
objects_list <- list(BPCells = Xen_Pros_bpcells, HDF5 = Xen_Pros_hdf5, HDF5_dense = Xen_Pros_hdf5_dense, Zarr = Xen_Pros_zarr)

# normalize
normalize.time <- list()
for(nm in names(objects_list)){
  start <- proc.time()
  objects_list[[nm]] <- normalizeData(objects_list[[nm]], sizefactor = 1000)
  normalize.time[[nm]] <- (proc.time() - start)[3]
}

# filter
for(nm in names(objects_list)){
  start <- proc.time()
  obj <- objects_list[[nm]]
  obj <- subset(obj, spatialpoints = vrSpatialPoints(obj)[obj$Count > 5])
  objects_list[[nm]] <- obj
}

# feature selection
feature.time <- list()
for(nm in names(objects_list)){
  print(nm)
  start <- proc.time()
  objects_list[[nm]] <- getFeatures(objects_list[[nm]], n = 2000)
  feature.time[[nm]] <- (proc.time() - start)[3]
}

# pca
pca.time <- list()
for(nm in names(objects_list)){
  print(nm)
  start <- proc.time()
  selected.features <- getVariableFeatures(objects_list[[nm]], n = 2000)
  objects_list[[nm]] <- getPCA(objects_list[[nm]], features = selected.features, dims = 30)
  pca.time[[nm]] <- (proc.time() - start)[3]
}

###
## results ####
###

results <- data.frame(savedisk = unlist(savedisk.time), 
                      normalize = unlist(normalize.time),
                      featureselection = unlist(feature.time),
                      pca = unlist(pca.time))
write.csv(results, file = "results/Xenium5000.csv", quote = FALSE, row.names = TRUE)

###
## Plot ####
###

results <- read.csv("results/Xenium5000.csv")
results_vis <- reshape2::melt(results, id.vars = "X")
ggplot(results_vis, mapping = aes(x = X, y = value, fill = X)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~variable, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  xlab("Method") + ylab("seconds")
ggsave(filename = "results/Xenium5000.pdf", plot = last_plot(), device = "pdf", width = 10, height = 8)