print("Loading packages")
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("ggplot2")) 
suppressPackageStartupMessages(library("dplyr"))

pcs <- 10
res <- 0.1
print(paste("using", pcs, "PCs"))
print(paste("using", res, "clustering resolution"))

print("Reading in integrated Seurat object")
integrated_samples <- readRDS(file="../../../seurat_objects/all_samples_integrated_2000features.RDS")

print("Running UMAP")
DefaultAssay(integrated_samples) <- "integrated"
integrated_samples <- RunUMAP(integrated_samples, reduction = "pca", dims = 1:pcs)
print("UMAP embedding completed!")

print("Clustering")
integrated_samples <- FindNeighbors(integrated_samples, reduction = "pca", dims = 1:pcs)
integrated_samples <- FindClusters(integrated_samples, resolution = res)
print("Clustering finished!")

print("Saving the reduced and clustered dataset")
integrated_fname <- paste0("../../../seurat_objects/all_samples_integrated_clustered_2000features_",pcs,"pcs_",res,"res.RDS")
saveRDS(integrated_samples, file=integrated_fname)

print(integrated_samples)
print("All done!")

