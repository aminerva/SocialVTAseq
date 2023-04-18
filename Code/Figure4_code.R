# Script for creating the plots in Figure 4

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
suppressPackageStartupMessages(library("NatParksPalettes"))         
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("reldist"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("tidyr"))

set.seed(305)

pcs <- 10 
res <- 0.3 

hungry_sated_pal <- natparks.pals("Glacier", 5)
hungry_sated_pal <- c(hungry_sated_pal[2], # sated
            		  "gray", # control
            		  hungry_sated_pal[4]) # hungry

da_pal <- c("#B2C247",
            "#1B876F",
            "#9EDF94",
            "#1B2E1A",
            "#778071",
            "#61C1A4",
            "#B6DF47",
            "#366C0B",
            "#A4FF67")

print("Reading in csvs with previously calculated DEGs")
da_degs <- read.csv(file="../../../analysis/da_degs.csv")
gaba_degs <- read.csv(file="../../../analysis/gaba_degs.csv")
glut_degs <- read.csv(file="../../../analysis/glut_degs.csv")
all_hungry_sated_degs <- read.csv(file="../../../analysis/all_hungry_sated_degs.csv")

######## panel B
print("Making volcano plot for each celltype")
cols <- c("down"=hungry_sated_pal[1],"no"=hungry_sated_pal[2],"up"=hungry_sated_pal[3])
options(repr.plot.width = 15, repr.plot.height = 5) 
a <- ggplot(da_degs, aes(x=avg_log2FC, y=-log10(p_val_adj), color=diffexpressed)) +
        geom_point(size = 0.6) + xlim(-1,1.2) + 
        scale_color_manual(values=cols) + 
        ggtitle("dopamine", "adjusted p<0.01 and log2FC>0.38") +
        labs(x="higher in sated <--   log2FC   --> higher in hungry", y="-10log(adj.p)") +
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
        geom_point(size = 0.6) + xlim(-1,1.2) + 
        scale_color_manual(values=cols) +
        ggtitle("gaba", "adjusted p<0.01 and log2FC>0.38") +
        labs(x="higher in sated <--   log2FC   --> higher in hungry", y="-10log(adj.p)") +
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
        geom_point(size = 0.6) + xlim(-1,1.2) + 
        scale_color_manual(values=cols) + 
        ggtitle("glut", "adjusted p<0.01 and log2FC>0.38") +
        labs(x="higher in sated <--   log2FC   --> higher in hungry", y="-10log(adj.p)") +
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
ggsave(file="../../../paper_figures/4B.pdf",
       width=15, height=5,dpi=320)

print("Making bar plot with DEG counts across celltypes")
cols <- c("dopamine"="#7BD1A9", "gaba"="#6BAA8D", "glut"="#456F5C")
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

######## panel E
print("Reading in DA clustered seurat object")
all_da <- readRDS(file=paste0("../../../seurat_objects/all_da_integrated_clustered_2000features_",pcs,"pcs_",res,"res.RDS"))

print("Plotting UMAP")
plotname2 <- paste0("../../../plots/all_da_integrated_clustered_2000features_",pcs,"pcs_",res,"res_umap.png")
DimPlot(all_da, reduction="umap", label=FALSE, pt.size=0.1, order=FALSE, repel=TRUE) + 
    scale_colour_manual(values = da_pal) + NoAxes() + NoLegend()
ggsave("../../../paper_figures/4E.pdf", width=4, height=3,dpi=320)

######## panel F

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
write.csv(da_top5_unique, file="../../../analysis/all_da_top5_unique_markers.csv")
DoHeatmap(all_da, features=c(da_top5_unique[['gene']],food_social_genes), slot="data",
          angle=0, group.colors=da_pal) +
    scale_fill_gradient2(low = '#2166ac', mid = "white", high = "#3E3E3E",
                         midpoint = 0, guide = "colourbar", aesthetics="fill") 
ggsave("../../../paper_figures/4F.pdf", width=10, height=3,dpi=320)

######## panel G
gene1 <- XX
print(paste("Making feature plot for", gene1, "in dopamine neurons")
if (gene1 %in% food_genes) {
        fp <- FeaturePlot(all_da, features=gene1, slot="data", 
                     order=TRUE, cols=c("#e8e8e8","#855e41"), 
                         pt.size=0.5) + labs(title=gene1) + NoAxes()
    } else if (food_social_genes[gene] %in% social_genes) {
        fp <- FeaturePlot(all_da, features=gene1, slot="data", 
                     order=TRUE, cols=c("#e8e8e8","#8BB1DD"), 
                         pt.size=0.5) + labs(title=gene1) + NoAxes()
    } 
savename1 <- paste0("../../../plots/",gene1,"_da_neurons_featureplot.pdf")
    ggsave(savename1, width=4, height=3,dpi=320)

######## panel H
gene2 <- XX
print(paste("Making feature plot for", gene2, "in dopamine neurons")
if (gene2 %in% food_genes) {
        fp <- FeaturePlot(all_da, features=gene2, slot="data", 
                     order=TRUE, cols=c("#e8e8e8","#855e41"), 
                         pt.size=0.5) + labs(title=gene2) + NoAxes()
} else if (gene2 %in% social_genes) {
    fp <- FeaturePlot(all_da, features=gene2, slot="data", 
                 order=TRUE, cols=c("#e8e8e8","#8BB1DD"), 
                     pt.size=0.5) + labs(title=gene2) + NoAxes()
} 
savename2 <- paste0("../../../plots/",gene2,"_da_neurons_featureplot.pdf")
    ggsave(savename2, width=4, height=3,dpi=320)
















