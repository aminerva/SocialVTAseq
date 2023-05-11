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

hungry_sated_pal <- natparks.pals("Glacier", 5)
hungry_sated_pal <- c(hungry_sated_pal[2], # sated
                      "gray", # control
                      hungry_sated_pal[4]) # hungry

set.seed(305)
print("Reading in the seurat objects")

hungry_sated_da <- readRDS(file="../../../seurat_objects/hungry_sated_da.RDS") 
hungry_sated_gaba <- readRDS(file="../../../seurat_objects/hungry_sated_gaba.RDS") 
hungry_sated_glut <- readRDS(file="../../../seurat_objects/hungry_sated_glut.RDS") 

Idents(hungry_sated_da) <- "condition"
Idents(hungry_sated_gaba) <- "condition"
Idents(hungry_sated_glut) <- "condition"

print("Finding calculating expression across conditions")
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

print("Assigning DEG labels based on logFC and pvalue")
da_degs$diffexpressed <- "no"
da_degs$diffexpressed[da_degs$avg_log2FC > 0.38 & da_degs$p_val_adj < 0.05] <- "up"
da_degs$diffexpressed[da_degs$avg_log2FC < -0.38 & da_degs$p_val_adj < 0.05] <- "down"

gaba_degs$diffexpressed <- "no"
gaba_degs$diffexpressed[gaba_degs$avg_log2FC > 0.38 & gaba_degs$p_val_adj < 0.05] <- "up"
gaba_degs$diffexpressed[gaba_degs$avg_log2FC < -0.38 & gaba_degs$p_val_adj < 0.05] <- "down"

glut_degs$diffexpressed <- "no"
glut_degs$diffexpressed[glut_degs$avg_log2FC > 0.38 & glut_degs$p_val_adj < 0.05] <- "up"
glut_degs$diffexpressed[glut_degs$avg_log2FC < -0.38 & glut_degs$p_val_adj < 0.05] <- "down"

all_hungry_sated_degs <- rbind(da_degs,gaba_degs,glut_degs)

print("Saving results to csv")
write.csv(da_degs, file="../../../analysis/da_degs.csv")
write.csv(gaba_degs, file="../../../analysis/gaba_degs.csv")
write.csv(glut_degs, file="../../../analysis/glut_degs.csv")
write.csv(all_hungry_sated_degs, file="../../../analysis/all_hungry_sated_degs.csv")

print("Making volcano plots")
cols <- c("down"=hungry_sated_pal[3],"no"=hungry_sated_pal[2],"up"=hungry_sated_pal[1])
options(repr.plot.width = 15, repr.plot.height = 5) 
a <- ggplot(da_degs, aes(x=avg_log2FC, y=-log10(p_val_adj), color=diffexpressed)) +
        geom_point(size = 0.1) + xlim(-1,1.2) + 
        scale_color_manual(values=cols) + 
        ggtitle("dopamine", "adjusted p<0.01 and log2FC>0.38") +
        labs(x="higher in hungry <--   log2FC   --> higher in sated", y="-10log(adj.p)") +
        geom_label_repel(label=da_degs$label, point.padding=unit(0.4,"lines"), 
                         nudge_x=0, nudge_y=0.5, label.size=0, fill=NA, 
                         min.segment.length=unit(0.1,"lines"), max.overlaps=50, force=2) + 
        NoLegend() + scale_y_continuous(expand=c(0,0), limits=c(-5,200)) +
        theme(axis.text.x = element_text(size=10, color="black"),
              axis.text.y = element_text(size=10, color="black"),
              axis.title.x = element_text(size=12, color="black"),
              axis.title.y = element_text(size=12, color="black"),
              plot.title = element_text(size=12))  
b <- ggplot(gaba_degs, aes(x=avg_log2FC, y=-log10(p_val_adj), color=diffexpressed)) +
        geom_point(size = 0.1) + xlim(-1,1.2) + 
        scale_color_manual(values=cols) +
        ggtitle("gaba", "adjusted p<0.01 and log2FC>0.38") +
        labs(x="higher in hungry <--   log2FC   --> higher in sated", y="-10log(adj.p)") +
        geom_label_repel(label=gaba_degs$label, point.padding=unit(0.4,"lines"), 
                         nudge_x=0, nudge_y=0.5, label.size=0, fill=NA, 
                         min.segment.length=unit(0.1,"lines"), max.overlaps=50, force=2)+ 
        NoLegend() + scale_y_continuous(expand=c(0,0), limits=c(-5,200)) +
        theme(axis.text.x = element_text(size=10, color="black"),
              axis.text.y = element_text(size=10, color="black"),
              axis.title.x = element_text(size=12, color="black"),
              axis.title.y = element_text(size=12, color="black"),
              plot.title = element_text(size=12))  
c <- ggplot(glut_degs, aes(x=avg_log2FC, y=-log10(p_val_adj), color=diffexpressed)) +
        geom_point(size = 0.1) + xlim(-1,1.2) + 
        scale_color_manual(values=cols) + 
        ggtitle("glut", "adjusted p<0.01 and log2FC>0.38") +
        labs(x="higher in hungry <--   log2FC   --> higher in sated", y="-10log(adj.p)") +
        geom_label_repel(label=glut_degs$label, point.padding=unit(0.4,"lines"), 
                         nudge_x=0, nudge_y=0.5, label.size=0, fill=NA, 
                         min.segment.length=unit(0.1,"lines"), max.overlaps=50, force=2) + 
        NoLegend() + scale_y_continuous(expand=c(0,0), limits=c(-5,200)) +
        theme(axis.text.x = element_text(size=10, color="black"),
              axis.text.y = element_text(size=10, color="black"),
              axis.title.x = element_text(size=12, color="black"),
              axis.title.y = element_text(size=12, color="black"),
              plot.title = element_text(size=12)) 
a | b | c
ggsave(file="../../../plots/hungry_sated_deg_volcanos.pdf",
       width=15, height=5,dpi=320)

print("Making bar plot with DEG counts across celltypes")
#cols <- c("dopamine"="#7BD1A9", "gaba"="#6BAA8D", "glut"="#456F5C")
cols <- c("dopamine"="#65A286", "gaba"="#965cb5", "glut"="#e8d568")
options(repr.plot.width = 3.5, repr.plot.height = 4) 
all_hungry_sated_degs <- all_hungry_sated_degs[all_hungry_sated_degs$p_val_adj<0.05 & abs(all_hungry_sated_degs$avg_log2FC)>0.38,]
ggplot(na.omit(all_hungry_sated_degs), aes(x=celltype, fill=celltype)) + geom_bar() + 
    scale_fill_manual(values=cols) +
    theme_classic() + labs(title="number of DEGs") +
    theme(axis.text.x = element_text(size=15, color="black"),
          axis.text.y = element_text(size=15, color="black"),
          axis.title.x = element_text(size=18, color="black"),
          axis.title.y = element_text(size=18, color="black"),
          plot.title = element_text(size=18, hjust=0.5)) + NoLegend() +
    scale_y_continuous(expand=c(0,0)) 
ggsave(file="../../../plots/hungry_sated_deg_counts.pdf",
       width=3.5, height=4,dpi=320)





