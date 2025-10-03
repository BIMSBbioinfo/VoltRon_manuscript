library(VoltRon)

###
# Xenium 313 panel ####
###

# import
Xen_R1 <- importXenium("../../../test/mainWorkflows/data/VisiumXenium/10X_Xenium_Visium/Xenium_R1/outs/", sample_name = "XeniumR1",
                       resolution_level = 3, overwrite_resolution = TRUE, import_molecules = FALSE)

###
## save to disk ####
###

# time
savedisk.time <- list()

# save with BPCells
start <- proc.time()
Xen_R1_bpcells <- saveVoltRon(Xen_R1,
                           output = "data/Xen_R1_bpcells/",
                           format = "HDF5VoltRon",
                           replace = TRUE,
                           verbose = FALSE)
savedisk.time[["BPCells"]] <- (proc.time() - start)[3]

# save with HDF5 DelayedArray
start <- proc.time()
Xen_R1_hdf5 <- saveVoltRon(Xen_R1,
                           output = "data/Xen_R1_hdf5/",
                           format = "HDF5VoltRon",
                           replace = TRUE,
                           verbose = FALSE, 
                           as.sparse = TRUE,
                           feature.vs.obs.engine = "DelayedArray")
savedisk.time[["HDF5"]] <- (proc.time() - start)[3]

# save with HDF5 DelayedArray
start <- proc.time()
Xen_R1_hdf5_dense <- saveVoltRon(Xen_R1,
                           output = "data/Xen_R1_hdf5_dense/",
                           format = "HDF5VoltRon",
                           replace = TRUE,
                           verbose = FALSE, 
                           as.sparse = FALSE,
                           feature.vs.obs.engine = "DelayedArray")
savedisk.time[["HDF5_dense"]] <- (proc.time() - start)[3]

# save with Zarr
start <- proc.time()
Xen_R1_zarr <- saveVoltRon(Xen_R1,
                           output = "data/Xen_R1_zarr/",
                           format = "ZarrVoltRon",
                           replace = TRUE,
                           verbose = FALSE)
savedisk.time[["Zarr"]] <- (proc.time() - start)[3]


###
## workflow ####
###

# load
Xen_R1_bpcells <- loadVoltRon(dir = "data/Xen_R1_bpcells/")
Xen_R1_hdf5 <- loadVoltRon(dir = "data/Xen_R1_hdf5/")
Xen_R1_hdf5_dense <- loadVoltRon(dir = "data/Xen_R1_hdf5_dense/")
Xen_R1_zarr <- loadVoltRon(dir = "data/Xen_R1_zarr/")
objects_list <- list(BPCells = Xen_R1_bpcells, HDF5 = Xen_R1_hdf5, HDF5_dense = Xen_R1_hdf5_dense, Zarr = Xen_R1_zarr)

# normalize
normalize.time <- list()
for(nm in names(objects_list)){
  start <- proc.time()
  objects_list[[nm]] <- normalizeData(objects_list[[nm]], sizefactor = 1000)
  normalize.time[[nm]] <- (proc.time() - start)[3]
}

# filter
# for(nm in names(objects_list)){
#   obj <- objects_list[[nm]]
#   obj <- subset(obj, spatialpoints = vrSpatialPoints(obj)[obj$Count > 5])
#   objects_list[[nm]] <- obj
# }

# temp
BPCells::write
normdata <- objects_list$BPCells[["Assay1"]]@data$RNA_norm
svd <- BPCells::svds(normdata, k=30, threads = 2L)

# pca
pca.time <- list()
for(nm in names(objects_list)){
  print(nm)
  start <- proc.time()
  n.workers <- if(nm == "BPCells") 1L else 3L
  objects_list[[nm]] <- getPCA(objects_list[[nm]], dims = 30, n.workers = n.workers)
  pca.time[[nm]] <- (proc.time() - start)[3]
}

# compare Rspectra, BiocSingular and BPCells
dims <- 30
pca_new.time <- list()
# bpcells time 
start <- proc.time()
normdata <- objects_list$BPCells[["Assay1"]]@data$RNA_norm
svd <- BPCells::svds(normdata, k=dims)
pr.data <- BPCells::multiply_cols(svd$v, svd$d)
pca_new.time[["BPCells"]] <- (proc.time() - start)[3]
# biocsingular
start <- proc.time()
normdata <- t(objects_list$HDF5[["Assay1"]]@data$RNA_norm)
svd <- BiocSingular::runSVD(normdata, k = dims)
x <- sweep(svd$u, 2, svd$d, "*")
pca_new.time[["BiocSingular"]] <- (proc.time() - start)[3]
# RSpectra
start <- proc.time()
mat <- objects_list$HDF5[["Assay1"]]@data$RNA_norm
row_means <- DelayedMatrixStats::rowMeans2(mat)
Ax <- function(x, args)
  (as.numeric(mat %*% x) - row_means * sum(x))
Atx <- function(x, args)
  (as.numeric(x %*% mat) - as.vector(row_means %*% x))
svd <- RSpectra::svds(Ax, Atrans=Atx, k=dims, dim=dim(mat))
# svd <- RSpectra::svds(mat, k=dims, dim=dim(mat))
x <- sweep(svd$u, 2, svd$d, "*")
pca_new.time[["RSpectra"]] <- (proc.time() - start)[3]

###
## results ####
###

results <- data.frame(savedisk = unlist(savedisk.time), 
                      normalize = unlist(normalize.time),
                      #featureselection = unlist(feature.time), 
                      pca = unlist(pca.time))
write.csv(results, file = "results/Xenium313.csv", quote = FALSE, row.names = TRUE)

# BPCells error
library(BPCells)
m <- matrix(data = seq_len(100*100), nrow=100) |>
  as("IterableMatrix")
file <- tempfile(fileext = ".h5")
m <- write_matrix_hdf5(m, path = file, group = "name")
# file <- tempdir()
# m <- write_matrix_dir(m, dir = file)
svd <- BPCells::svds(m, k=30, threads=2L)
