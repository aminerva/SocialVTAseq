
print("Loading packages")
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("ggplot2")) 
suppressPackageStartupMessages(library("dplyr"))

print("Reading in Seurat objects")
sample_names <- c("MCVTA1", "MCVTA2", "MCVTA3", "FCVTA1", "FCVTA2", "FCVTA3", 
                  "hungryM", "hungryF", "satedM", "satedF")
for (i in 1:length(sample_names)) {
    filename <- paste0("../../../seurat_objects/", sample_names[i], "_preprocessed.RDS")
    obj <- readRDS(filename)
    assign(sample_names[i], obj)
}
samples <- c(MCVTA1, MCVTA2, MCVTA3, FCVTA1, FCVTA2, FCVTA3, hungryM, hungryF, satedM, satedF)

nFeats <- 2000
pcs <- 10
res <- 0.1

print("Normalzing, finding variable features, scaling, and running PCA on each sample")

# Alternative workflow for integrating large datasets consists of the following steps:
    # Create a list of Seurat objects to integrate
    # Perform normalization, feature selection, and scaling separately for each dataset
    # Run PCA on each object in the list
    # Integrate datasets, and proceed with joint analysis

print("Normalizing and finding variable features")
samples <- lapply(X = samples, FUN = function(x) {
    # The fx NormalizeData() takes the feature counts for each cell and divides them by 
    # the total counts for that cell. These feature counts are then multiplied by the 
    # scale.factor and natural-log transformed using log1p.
    # THE OUTPUT OF THIS IS EXPRESSION IN TRANSCRIPTS PER 10,000
    x <- NormalizeData(x, verbose = TRUE, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nFeats, verbose = TRUE)
})

print("Selecting features that are repeatedly variable across the datasets for integration")
int.features <- SelectIntegrationFeatures(object.list = samples, nfeatures = nFeats)
fname1 <- paste0('../../../seurat_objects/integration_features_',nFeats,'features.RDS')
saveRDS(int.features, file = fname1)
print("Found integration features")

print("Scaling and running PCA on each object separately")
samples <- lapply(X = samples, FUN = function(x) {
    x <- ScaleData(x, features = int.features, verbose = TRUE)
    x <- RunPCA(x, features = int.features, verbose = TRUE)
})

print("Saving each separate normalized and scaled Seurat object")
for (i in 1:length(samples)) {
    savename <- paste0("../../../seurat_objects/", sample_names[i], "_normalized.RDS")
    saveRDS(samples[[i]], savename)
}


print("Integrating all datasets")

print("Finding integration anchors")
int.anchors <- FindIntegrationAnchors(object.list=samples, reduction="rpca", 
                                      dims=1:50, anchor.features=nFeats)
fname2 <- paste0('../../../seurat_objects/integration_anchors_',nFeats,'features.RDS')
saveRDS(int.anchors, file=fname2)
print("Found integration anchors!")

print("Integrating datasets")
integrated_samples <- IntegrateData(anchorset = int.anchors, dims=1:50)
"Datasets successfully integrated!"

# Specify that downstream analysis is on integrated data 
DefaultAssay(integrated_samples) <- "integrated" 

print("Scaling and running PCA on the integrated dataset")
integrated_samples <- ScaleData(integrated_samples)
integrated_samples <- RunPCA(integrated_samples, npcs = 50)

nuclei_num <- length(Cells(integrated_samples))
print(paste0("The integrated dataset has ", nuclei_num, " total nuclei"))

print("Saving")
fname3 <- paste0('../../../seurat_objects/all_samples_integrated_',nFeats,'features.RDS')
saveRDS(integrated_samples, file=fname3)

ElbowPlot(integrated_samples, ndim=50)
plotname1 <- paste0("../../../plots/all_samples_integrated_",nFeats,"features_elbow.png")
ggsave(plotname1, width=10, height=6,dpi=320)

print("All done!")
print("Look at elbow plot for PCs to use for UMAP embedding and clsutering")

