print("Loading packages")
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("ggplot2")) 
suppressPackageStartupMessages(library("dplyr"))

print("Reading in parameters (dimensions (PCs) for UMAP embedding and resolution for clustering)")
args <- commandArgs(trailingOnly = TRUE)
pcs <- as.numeric(args[1])
res <- as.numeric(args[2])
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


print("Plotting UMAP and DotPlot")
DefaultAssay(integrated_samples) <- "RNA"
print("Making UMAP")
plotname2 <- paste0("../../../plots/all_samples_integrated_clustered_2000features_",pcs,"pcs_",res,"res_umap.png")
DimPlot(integrated_samples, reduction = "umap", label = FALSE, repel = TRUE)
ggsave(plotname2, width=6, height=6,dpi=320)


print("Making dotplot")
marker_features <- rev(c("Syt1", "Rbfox3", 
                         "Th", "Ddc", "Slc6a3",
                         "Slc17a6", "Gad1", "Gad2", 
                         "Aqp4", "Cx3cr1", "Pdgfra", "Mog", "Mobp", "Cldn5", "Vtn"))
plotname3 <- paste0("../../../plots/all_samples_integrated_clustered_2000features_",pcs,"pcs_",res,"res_dotplot.png")
DotPlot(integrated_samples, features = marker_features, dot.scale=10, col.min=0, 
        cols=c("#F9F9F9","darkorchid4")) + labs(y="Cluster", x="Gene") + coord_flip() +
        theme(plot.background=element_rect(fill="white", color="white"))
ggsave(plotname3, width=10, height=6,dpi=320)

print(integrated_samples)
print("All done!")

