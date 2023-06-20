print("Loading packages")
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("ggplot2")) 
suppressPackageStartupMessages(library("dplyr"))

pcs <- 10
res <- 0.3
print(paste("using", pcs, "PCs"))
print(paste("using", res, "clustering resolution"))

print("Reading in the neurons only Seurat object")
all_neurons <- readRDS(file="../../../seurat_objects/all_conditions_all_neurons.RDS")
DefaultAssay(all_neurons) <- "RNA"

print("Pulling out just dopamine neurons")
all_da <- subset(all_neurons, subset = Slc6a3 > 0 | Th > 0 ) 
all_da[["nucleusID"]] <- colnames(all_da)
saveRDS(all_da, file="../../../seurat_objects/all_da.RDS") 

print("Finding variable features, scaling, and running PCA")
all_da <- FindVariableFeatures(all_da, selection.method="vst", nfeatures=2000)
all_da <- ScaleData(all_da)
all_da <- RunPCA(all_da, npcs = 50)

print("Running UMAP")
all_da <- RunUMAP(all_da, reduction = "pca", dims = 1:pcs, n.neighbors=70L, 
                   min.dist=0.3, spread=0.7)
all_da <- FindNeighbors(all_da, reduction = "pca", dims = 1:pcs, k.param=70)
print("UMAP embedding completed!")

print("Clustering")
all_da <- FindClusters(all_da, resolution = res)
print("Clustering finished!")

print("Saving the reduced and clustered dataset")
all_da_fname <- paste0("../../../seurat_objects/all_da_integrated_clustered_2000features_",pcs,"pcs_",res,"res.RDS")
saveRDS(all_da, file=all_da_fname)


print("Finding marker genes for each DA subcluster")
Idents(all_da) <- all_da[["seurat_clusters"]]
all_da_markers <- FindAllMarkers(all_da, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(all_da_markers, file=paste0("../../../analysis/all_da_",pcs,"pcs_",res,"res_markers.csv"))

print("All done!")

