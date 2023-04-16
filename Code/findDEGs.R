print("Loading packages")
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
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("tidyr"))

set.seed(305)
print("Reading in the seurat objects")

hungry_sated_da <- readRDS(file="../../../seurat_objects/hungry_sated_da.RDS") 
hungry_sated_gaba <- readRDS(file="../../../seurat_objects/hungry_sated_gaba.RDS") 
hungry_sated_glut <- readRDS(file="../../../seurat_objects/hungry_sated_glut.RDS") 

Idents(hungry_sated_da) <- "condition"
Idents(hungry_sated_gaba) <- "condition"
Idents(hungry_sated_glut) <- "condition"

da_degs <- FindMarkers(hungry_sated_da, ident.1="sated", ident.2="hungry", test.use="DESeq2",
                              logfc.threshold=0)
gaba_degs <- FindMarkers(hungry_sated_gaba, ident.1="sated", ident.2="hungry", test.use="DESeq2",
                                logfc.threshold=0)
glut_degs <- FindMarkers(hungry_sated_glut, ident.1="sated", ident.2="hungry", test.use="DESeq2",
                                logfc.threshold=0)

da_degs$label <- ifelse(abs(da_degs$avg_log2FC)>0.38 & da_degs$p_val_adj<0.05, rownames(da_degs), NA)
gaba_degs$label <- ifelse(abs(gaba_degs$avg_log2FC)>0.38 & gaba_degs$p_val_adj<0.05, rownames(gaba_degs), NA)
glut_degs$label <- ifelse(abs(glut_degs$avg_log2FC)>0.38 & glut_degs$p_val_adj<0.05, rownames(glut_degs), NA)

da_degs$celltype <- rep("dopamine",nrow(da_degs))
gaba_degs$celltype <- rep("gaba",nrow(gaba_degs))
glut_degs$celltype <- rep("glut",nrow(glut_degs))

all_hungry_sated_degs <- rbind(da_degs,gaba_degs,glut_degs)

write.csv(da_degs, file="../../../analysis/da_degs.csv")
write.csv(gaba_degs, file="../../../analysis/gaba_degs.csv")
write.csv(glut_degs, file="../../../analysis/glut_degs.csv")
write.csv(all_hungry_sated_degs, file="../../../analysis/all_hungry_sated_degs.csv")





