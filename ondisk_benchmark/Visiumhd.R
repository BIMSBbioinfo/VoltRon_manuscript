library(VoltRon)
library(reshape2)
library(ggplot2)

###
# Visium HD ####
###

# import
VisiumHD_exp <- importVisiumHD("../../../../data/visiumhd/10x_human_kidney/", bin.size = "8", 
                               sample_name = "VisiumHD_R1", 
                               resolution_level = "hires", inTissue = TRUE)

###
## save to disk ####
###

# time
savedisk.time <- list()

# save with BPCells
print("BPCells")
start <- proc.time()
VisiumHD_exp_bpcells <- saveVoltRon(VisiumHD_exp,
                                    output = "../data/ondisk_benchmark/VisiumHD_exp_bpcells/",
                                    format = "HDF5VoltRon",
                                    replace = TRUE,
                                    verbose = FALSE)
savedisk.time[["BPCells"]] <- (proc.time() - start)[3]
print(savedisk.time[["BPCells"]])

# save with HDF5 DelayedArray
print("HDF5")
start <- proc.time()
VisiumHD_exp_hdf5 <- saveVoltRon(VisiumHD_exp,
                                 output = "../data/ondisk_benchmark/VisiumHD_exp_hdf5/",
                                 format = "HDF5VoltRon",
                                 replace = TRUE,
                                 verbose = FALSE, 
                                 as.sparse = TRUE,
                                 feature.vs.obs.engine = "DelayedArray")
savedisk.time[["HDF5"]] <- (proc.time() - start)[3]
print(savedisk.time[["HDF5"]])

# save with HDF5 DelayedArray
print("HDF5_dense")
start <- proc.time()
VisiumHD_exp_hdf5_dense <- saveVoltRon(VisiumHD_exp,
                                 output = "../data/ondisk_benchmark/VisiumHD_exp_hdf5_dense/",
                                 format = "HDF5VoltRon",
                                 replace = TRUE,
                                 verbose = FALSE, 
                                 as.sparse = FALSE,
                                 feature.vs.obs.engine = "DelayedArray")
savedisk.time[["HDF5_dense"]] <- (proc.time() - start)[3]
print(savedisk.time[["HDF5_dense"]])

# save with Zarr
print("Zarr")
start <- proc.time()
VisiumHD_exp_zarr <- saveVoltRon(VisiumHD_exp,
                           output = "../data/ondisk_benchmark/VisiumHD_exp_zarr/",
                           format = "ZarrVoltRon",
                           replace = TRUE,
                           verbose = FALSE)
savedisk.time[["Zarr"]] <- (proc.time() - start)[3]
print(savedisk.time[["Zarr"]])

###
## workflow ####
###

# load
VisiumHD_exp_bpcells <- loadVoltRon(dir = "../data/ondisk_benchmark/VisiumHD_exp_bpcells/")
VisiumHD_exp_hdf5 <- loadVoltRon(dir = "../data/ondisk_benchmark/VisiumHD_exp_hdf5/")
VisiumHD_exp_hdf5_dense <- loadVoltRon(dir = "../data/ondisk_benchmark/VisiumHD_exp_hdf5_dense/")
VisiumHD_exp_zarr <- loadVoltRon(dir = "../data/ondisk_benchmark/VisiumHD_exp_zarr/")
objects_list <- list(BPCells = VisiumHD_exp_bpcells, HDF5 = VisiumHD_exp_hdf5, HDF5_dense = VisiumHD_exp_hdf5_dense, Zarr = VisiumHD_exp_zarr)

# normalize
normalize.time <- list()
message("normalization")
for(nm in names(objects_list)){
  print(nm)
  start <- proc.time()
  objects_list[[nm]] <- normalizeData(objects_list[[nm]], sizefactor = 1000)
  normalize.time[[nm]] <- (proc.time() - start)[3]
  print(normalize.time[[nm]])
}
print(unlist(normalize.time))

# filter
for(nm in names(objects_list)){
  start <- proc.time()
  obj <- objects_list[[nm]]
  obj <- subset(obj, spatialpoints = vrSpatialPoints(obj)[obj$Count > 5])
  objects_list[[nm]] <- obj
}

# feature selection
feature.time <- list()
message("feature selection")
for(nm in names(objects_list)){
  print(nm)
  start <- proc.time()
  objects_list[[nm]] <- getFeatures(objects_list[[nm]], n = 3000)
  feature.time[[nm]] <- (proc.time() - start)[3]
  print(feature.time[[nm]])
}
print(unlist(feature.time))

# pca
pca.time <- list()
message("pca reduction")
for(nm in names(objects_list)){
  print(nm)
  start <- proc.time()
  selected.features <- getVariableFeatures(objects_list[[nm]], n = 3000)
  objects_list[[nm]] <- getPCA(objects_list[[nm]], features = selected.features, dims = 30)
  pca.time[[nm]] <- (proc.time() - start)[3]
  print(pca.time[[nm]])
}
print(unlist(pca.time))

###
## results ####
###

results <- data.frame(savedisk = unlist(savedisk.time), 
                      normalize = unlist(normalize.time),
                      feature = unlist(feature.time),
                      pca = unlist(pca.time))
write.csv(results, file = "results/VisiumHD_exp.csv", quote = FALSE, row.names = TRUE)

###
## Plot ####
###

results <- read.csv("results/VisiumHD_exp.csv")
results_vis <- reshape2::melt(results, id.vars = "X")
ggplot(results_vis, mapping = aes(x = X, y = value, fill = X)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~variable, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  xlab("Method") + ylab("seconds")
ggsave(filename = "results/VisiumHD_exp.pdf", plot = last_plot(), device = "pdf", width = 10, height = 8)