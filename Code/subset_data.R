# Script for reading in the neuron seurat objects and subsetting them into celltype specific objects

suppressPackageStartupMessages(library("abind"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggbreak"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("grid"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("gtools"))
suppressPackageStartupMessages(library("likert"))
suppressPackageStartupMessages(library("Matrix"))
suppressPackageStartupMessages(library("metap"))          
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("reldist"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("sctransform"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("tidyr"))

set.seed(305)

# Read in the neuron-only seurat object from the control mice
ctrl_neurons <- readRDS(file="./CONTROL_food_social_vta/ctrl_neurons.RDS") 

# Split the neuron subtypes
da <- subset(ctrl_neurons, subset = Slc6a3 > 0 | Th > 0 ) 
gaba <- subset(ctrl_neurons, subset = Slc32a1>0 & Slc6a3==0 & Th==0 & Slc17a6==0 & Slc17a7==0 & Slc17a8==0 | Gad1>0 & Slc6a3==0 & Th==0 & Slc17a6==0 & Slc17a7==0 & Slc17a8==0 | Gad2>0 & Slc6a3==0 & Th==0 & Slc17a6==0 & Slc17a7==0 & Slc17a8==0)
glut <- subset(ctrl_neurons, subset = Slc17a6>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0 | Slc17a7>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0 | Slc17a8>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0) 

da[["nucleusID"]] <- colnames(da)
gaba[["nucleusID"]] <- colnames(gaba)
glut[["nucleusID"]] <- colnames(glut)

saveRDS(da, file="../social_vta/seurat_objects/ctrl_da_obj.RDS") 
saveRDS(gaba, file="../social_vta/seurat_objects/ctrl_gaba_obj.RDS") 
saveRDS(glut, file="../social_vta/seurat_objects/ctrl_glut_obj.RDS") 


# Read in data from hungry/sated experiment
# hungry_sated <- readRDS(file="../hungry_sated/integrated_samples_2000features_10pcs_0.05res.RDS")
# DefaultAssay(hungry_sated) <- "RNA"
# saveRDS(hungry_sated, file="../hungry_sated/hungry_sated.RDS")

# Pull out only neurons (clusters 1 and 3)
# hungry_sated_neurons <- subset(hungry_sated, idents = c("1", "3"))
# saveRDS(hungry_sated_neurons, file="../hungry_sated/hungry_sated_neurons.RDS") 
hungry_sated_neurons <- readRDS(file="../hungry_sated/hungry_sated_neurons.RDS") 

# Split the neuron subtypes
hungry_sated_da <- subset(hungry_sated_neurons, subset = Slc6a3 > 0 | Th > 0 ) 
hungry_sated_gaba <- subset(hungry_sated_neurons, subset = Slc32a1>0 & Slc6a3==0 & Th==0 & Slc17a6==0 & Slc17a7==0 & Slc17a8==0 | Gad1>0 & Slc6a3==0 & Th==0 & Slc17a6==0 & Slc17a7==0 & Slc17a8==0 | Gad2>0 & Slc6a3==0 & Th==0 & Slc17a6==0 & Slc17a7==0 & Slc17a8==0)
hungry_sated_glut <- subset(hungry_sated_neurons, subset = Slc17a6>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0 | Slc17a7>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0 | Slc17a8>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0) 

hungry_sated_da[["nucleusID"]] <- colnames(hungry_sated_da)
hungry_sated_gaba[["nucleusID"]] <- colnames(hungry_sated_gaba)
hungry_sated_glut[["nucleusID"]] <- colnames(hungry_sated_glut)

saveRDS(hungry_sated_da, file="../social_vta/seurat_objects/hungry_sated_da_obj.RDS") 
saveRDS(hungry_sated_gaba, file="../social_vta/seurat_objects/hungry_sated_gaba_obj.RDS") 
saveRDS(hungry_sated_glut, file="../social_vta/seurat_objects/hungry_sated_glut_obj.RDS") 

# Split by condition (hungry or sated)
hungry_da <- subset(hungry_sated_da, subset = condition == "hungry")
hungry_gaba <- subset(hungry_sated_gaba, subset = condition == "hungry")
hungry_glut <- subset(hungry_sated_glut, subset = condition == "hungry")
sated_da <- subset(hungry_sated_da, subset = condition == "sated")
sated_gaba <- subset(hungry_sated_gaba, subset = condition == "sated")
sated_glut <- subset(hungry_sated_glut, subset = condition == "sated")

saveRDS(hungry_da, file="../social_vta/seurat_objects/hungry_da_obj.RDS") 
saveRDS(hungry_gaba, file="../social_vta/seurat_objects/hungry_gaba_obj.RDS") 
saveRDS(hungry_glut, file="../social_vta/seurat_objects/hungry_glut_obj.RDS") 
saveRDS(sated_da, file="../social_vta/seurat_objects/sated_da_obj.RDS") 
saveRDS(sated_gaba, file="../social_vta/seurat_objects/sated_gaba_obj.RDS") 
saveRDS(sated_glut, file="../social_vta/seurat_objects/sated_glut_obj.RDS") 





