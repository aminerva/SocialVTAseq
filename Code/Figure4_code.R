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

hungry_sated_pal <- natparks.pals("Glacier", 5)
hungry_sated_pal <- c(hungry_sated_pal[2], # sated
            		  "gray", # control
            		  hungry_sated_pal[4]) # hungry

print("Reading in csvs with previously calculated DEGs")
da_degs <- read.csv(file="../../../analysis/da_degs.csv")
gaba_degs <- read.csv(file="../../../analysis/gaba_degs.csv")
glut_degs <- read.csv(file="../../../analysis/glut_degs.csv")
all_hungry_sated_degs <- read.csv(file="../../../analysis/all_hungry_sated_degs.csv")

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