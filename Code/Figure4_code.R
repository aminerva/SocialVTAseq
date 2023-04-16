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

da_degs$diffexpressed <- "no"
da_degs$diffexpressed[da_degs$avg_log2FC > 0.38 & da_degs$p_val_adj < 0.05] <- "up"
da_degs$diffexpressed[da_degs$avg_log2FC < 0.38 & da_degs$p_val_adj < 0.05] <- "down"


options(repr.plot.width = 15, repr.plot.height = 5) 
a <- ggplot(da_degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(size = 0.6, aes(color = diffexpressed)) + 
        scale_color_manual(values=hungry_sated_pal) + 
        scale_fill_manual(values=hungry_sated_pal) + xlim(-1,1.2) + 
        ggtitle("dopamine", "adjusted p<0.01 and log2FC>0.5") +
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

b <- ggplot(gaba_degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(size = 0.6, aes(color = p_val_adj < 0.05 & abs(avg_log2FC)> 0.38)) + 
        scale_color_manual(values=c("gray50","indianred")) + xlim(-1,1.2) + 
        ggtitle("gaba", "adjusted p<0.01 and log2FC>0.5") +
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

c <- ggplot(glut_degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(size = 0.6, aes(color = p_val_adj < 0.05 & abs(avg_log2FC)> 0.38)) + 
        scale_color_manual(values=c("gray50","indianred")) + xlim(-1,1.2) + 
        ggtitle("glut", "adjusted p<0.01 and log2FC>0.5") +
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

ggsave(file="../../../plots/hungry_sated_degs.pdf",
       width=15, height=5,dpi=320)