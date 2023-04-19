# Script for reading in the neuron seurat objects and subsetting them into celltype specific objects

print("Loading packages")
suppressPackageStartupMessages(library("abind"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Matrix"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("tidyr"))

set.seed(305)

print("Reading in the integrated object")
all_samples <- readRDS(file="../../../seurat_objects/all_samples_integrated_clustered_2000features_10pcs_0.05res.RDS") 
all_samples <- subset(all_samples, idents = c("0","1","2","3","4","5","6","7")) # Cluster 8 bad qual
DefaultAssay(all_samples) <- "RNA"

print("Pulling out neurons only")
all_neurons <- subset(all_samples, idents = c("1","3","6") )
DefaultAssay(all_neurons) <- "RNA"
saveRDS(all_neurons, file="../../../seurat_objects/all_conditions_all_neurons.RDS")

print("Creating separate objects for each condition (control, hungry, and sated)")
ctrl_neurons <- subset(all_neurons, subset = condition == "control")
hungry_sated_neurons <- subset(all_neurons, subset = condition %in% c("hungry","sated"))
saveRDS(ctrl_neurons, file="../../../seurat_objects/control_all_neurons.RDS")
saveRDS(hungry_sated_neurons, file="../../../seurat_objects/hungry_sated_all_neurons.RDS")

print("Creating separate objects for each celltype within each condition")
ctrl_da <- subset(ctrl_neurons, subset = Slc6a3 > 0 | Th > 0 ) 
ctrl_gaba <- subset(ctrl_neurons, subset = Slc32a1>0 & Slc6a3==0 & Th==0 & Slc17a6==0 & Slc17a7==0 & Slc17a8==0 | Gad1>0 & Slc6a3==0 & Th==0 & Slc17a6==0 & Slc17a7==0 & Slc17a8==0 | Gad2>0 & Slc6a3==0 & Th==0 & Slc17a6==0 & Slc17a7==0 & Slc17a8==0)
ctrl_glut <- subset(ctrl_neurons, subset = Slc17a6>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0 | Slc17a7>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0 | Slc17a8>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0) 
hungry_sated_da <- subset(hungry_sated_neurons, subset = Slc6a3 > 0 | Th > 0 ) 
hungry_sated_gaba <- subset(hungry_sated_neurons, subset = Slc32a1>0 & Slc6a3==0 & Th==0 & Slc17a6==0 & Slc17a7==0 & Slc17a8==0 | Gad1>0 & Slc6a3==0 & Th==0 & Slc17a6==0 & Slc17a7==0 & Slc17a8==0 | Gad2>0 & Slc6a3==0 & Th==0 & Slc17a6==0 & Slc17a7==0 & Slc17a8==0)
hungry_sated_glut <- subset(hungry_sated_neurons, subset = Slc17a6>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0 | Slc17a7>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0 | Slc17a8>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0) 


ctrl_da[["nucleusID"]] <- colnames(ctrl_da)
ctrl_gaba[["nucleusID"]] <- colnames(ctrl_gaba)
ctrl_glut[["nucleusID"]] <- colnames(ctrl_glut)
hungry_sated_da[["nucleusID"]] <- colnames(hungry_sated_da)
hungry_sated_gaba[["nucleusID"]] <- colnames(hungry_sated_gaba)
hungry_sated_glut[["nucleusID"]] <- colnames(hungry_sated_glut)

saveRDS(ctrl_da, file="../../../seurat_objects/ctrl_da.RDS") 
saveRDS(ctrl_gaba, file="../../../seurat_objects/ctrl_gaba.RDS") 
saveRDS(ctrl_glut, file="../../../seurat_objects/ctrl_glut.RDS") 
saveRDS(hungry_sated_da, file="../../../seurat_objects/hungry_sated_da.RDS") 
saveRDS(hungry_sated_gaba, file="../../../seurat_objects/hungry_sated_gaba.RDS") 
saveRDS(hungry_sated_glut, file="../../../seurat_objects/hungry_sated_glut.RDS") 

print("All done!")




