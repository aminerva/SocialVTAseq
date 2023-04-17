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


print("Plotting UMAP")
plotname2 <- paste0("../../../plots/all_da_integrated_clustered_2000features_",pcs,"pcs_",res,"res_umap.png")
DimPlot(all_da, reduction="umap", label=TRUE, pt.size=0.5, order=TRUE, repel=TRUE) + NoAxes() + NoLegend() 
    # scale_colour_manual(values = da_pal) #breaks = names(cell_counts), labels = cluster_labs, 
ggsave(plotname2, width=6, height=3,dpi=320)

print("Finding marker genes for each DA subcluster")
Idents(all_da) <- all_da[["seurat_clusters"]]
all_da_markers <- FindAllMarkers(all_da, only.pos = TRUE, logfc.threshold = 0.25)
all_da_top20 <- all_da_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
all_da_top10 <- all_da_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
all_da_top5 <- all_da_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(all_da_markers, 
          file=paste0("../../../analysis/all_da_",pcs,"pcs_",res,"res_markers.csv"))
write.csv(all_da_top20, 
          file=paste0("../../../analysis/all_da_",pcs,"pcs_",res,"res_top20_markers.csv"))
write.csv(all_da_top10, 
          file=paste0("../../../analysis/all_da_",pcs,"pcs_",res,"res_top10_markers.csv"))
write.csv(all_da_top5, 
          file=paste0("../../../analysis/all_da_",pcs,"pcs_",res,"res_top5_markers.csv"))

# List of genes relevant to feeding or social hormone receptors
food_receptors <- c("Cckar", "Cckbr","Ghsr", "Gipr", "Gpr171", "Hcrtr1", "Hcrtr2", "Insr",  "Lepr", 
                    "Mc1r", "Mc3r", "Mc4r", "Npy1r", "Pparg", "Slc16a2", "Thra", "Thrb", "Trhr", "Tshr") 
food_enzymes <- c("Pcsk1", "Pcsk2", "Cpe")
social_receptors <- c("Ar", "Avpr1a", "Esr1", "Esr2",  "Kiss1r", "Oxtr", "Pgr", "Prlr")  
social_enzymes <- c("Cyp19a1", "Cyp11a1", "Hsd3b2", "Srd5a1", "Srd5a2", "Srd5a3") 
food_genes <- c(food_receptors, food_enzymes)
social_genes <- c(social_receptors, social_enzymes)
food_social_genes <- c(food_genes, social_genes)

print("Plotting heatmap of DA subcluster marker genes")
da_markers <- read.csv(paste0("../../../analysis/all_da_",pcs,"pcs_",res,"res_markers.csv"))
da_top5_unique <- da_markers %>% ungroup() %>% distinct(gene, .keep_all = TRUE) %>% 
                           group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(da_top5_unique, file=paste0("../../../analysis/all_da_",pcs,"pcs_",res,"res_top5_unique_markers.csv"))

DoHeatmap(all_da, features=c(da_top5_unique[['gene']],food_social_genes), slot="data", angle=0) +
          scale_fill_gradient2(low = '#2166ac', mid = "white", high = "#3E3E3E",
                               midpoint = 0, guide = "colourbar", aesthetics="fill") 

ggsave(paste0("../../../plots/DA_markers_heatmap_",pcs,"pcs_",res,"res.pdf"), width=12, height=4,dpi=320)")

print("All done!")

