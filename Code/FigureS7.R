# Script for creating the plots in Figure S7

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
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("RColorBrewer"))    
suppressPackageStartupMessages(library("reldist"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("tidyr"))

source("./save_heatmap.R")

set.seed(305)

cluster_pal <- c('#F49FA8', # astro coral
	             '#D493EA', # endo purple
	             '#F4CB57', # micro gold 
	             '#7BC0A0', # neuron green
	             '#7B9CC0', # oligo blue
                 '#93D4EA', # OPC light blue/turquoise
                 '#7BC07D') #8 green
sex_pal <- c("#514289", "#C4878C")
hungry_sated_pal <- natparks.pals("Glacier", 5)
hungry_sated_pal <- c(hungry_sated_pal[2], # sated
                      hungry_sated_pal[4]) # hungry

print("Reading in integrated Seurat object")
fname <- paste0("../../../seurat_objects/all_samples_integrated_clustered_2000features_10pcs_0.05res.RDS")
integrated_samples <- readRDS(file=fname)
integrated_samples <- subset(integrated_samples, idents = c("0","1","2","3","4","5","6","7")) # Cluster 8 bad quality
DefaultAssay(integrated_samples) <- "RNA"

print("Adding cell name column to metadata")
celltype2 <- c()
for (cell in 1:nrow(integrated_samples[["seurat_clusters"]])) {
    if ((integrated_samples[["seurat_clusters"]][cell,]) %in% c("1","3","6")) {
        type <- "neuron"
    } else if ((integrated_samples[["seurat_clusters"]][cell,]) == "2") {
            type <- "astro"
    } else if ((integrated_samples[["seurat_clusters"]][cell,]) == "5") {
            type <- "micro"
    } else if ((integrated_samples[["seurat_clusters"]][cell,]) == "4") {
        type <- "OPC"
    } else if ((integrated_samples[["seurat_clusters"]][cell,]) == "0") {
        type <- "oligo"
    } else if ((integrated_samples[["seurat_clusters"]][cell,]) == "7") {
        type <- "endo"
    }
    celltype2 <- c(celltype2, type)
}
integrated_samples[["celltype_names"]] <- celltype2

print("Making plots")

######## panel A
print("Making UMAP") 
DimPlot(integrated_samples, reduction = "umap", label = FALSE, repel = TRUE, raster = FALSE,
		pt.size=0.1, group.by="celltype_names", cols=c(cluster_pal)) + labs(title="") + NoAxes() + NoLegend() 
ggsave(file="../../../paper_figures/S7A.pdf", width=5, height=5,dpi=320)

######## panel B
print("Making dotplot") 
Idents(integrated_samples) <- integrated_samples$celltype_names
Idents(integrated_samples) <- factor(Idents(integrated_samples), 
									 levels=c("neuron","astro","micro","OPC","oligo","endo"))
levels(integrated_samples) <- c("neuron","astro","micro","OPC","oligo","endo")
marker_features <- rev(c("Syt1",
                         #"Th", "Slc6a3",
                         #"Slc17a6", "Gad2", 
                         "Aqp4", "Flt1", "Cx3cr1", "Pdgfra", "Mobp"))
DotPlot(integrated_samples, features = marker_features, dot.scale=10, #col.min=0, 
        cols=c("white","#7BC0A0")) + labs(y="Cluster", x="Gene") + coord_flip() +
        theme(plot.background=element_rect(fill="white", color="white"),
              plot.margin = unit(c(1,0,0,0), "cm"),
              axis.text.x = element_text(angle=30, size=15, color="black", hjust=1, vjust=1),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size=15, color="black")) +
        scale_y_discrete(limits=c("neuron","astro","endo","micro","OPC","oligo"))
ggsave(file="../../../paper_figures/S7B.pdf", width=5.5, height=4,dpi=320)

######## panel C
print("Making violin plots") 
UMIs <- ggplot(data=integrated_samples@meta.data, aes(x=celltype_names, y=nCount_RNA, fill=celltype_names)) + 
        geom_violin(show.legend=F, trim=F, lwd=0.8, scale="width") + NoLegend() +
        scale_fill_manual(values=cluster_pal) +
        geom_boxplot(width=0.07, fill="white", outlier.size=0, lwd=0.8) + 
        theme_classic() + ylab("# UMIs") + 
        theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size=15, color="black"),
            axis.title.y = element_text(size=15, color="black")) +
        scale_x_discrete(limits=c("neuron","astro","endo","micro","OPC","oligo"))
genes <- ggplot(data=integrated_samples@meta.data, aes(x=celltype_names, y=nFeature_RNA, fill=celltype_names)) + 
        geom_violin(show.legend=F, trim=F, lwd=0.8, scale="width") + NoLegend() +
        scale_fill_manual(values=cluster_pal) +
        geom_boxplot(width=0.07, fill="white", outlier.size=0, lwd=0.8) + 
        theme_classic() + ylab("# genes") + 
        theme(axis.text.x = element_text(size=15, color="black", angle=30, vjust=1, hjust=1),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size=15, color="black"),
            axis.title.y = element_text(size=15, color="black")) +
        scale_x_discrete(limits=c("neuron","astro","endo","micro","OPC","oligo"))
ggsave(UMIs / genes, file="../../../paper_figures/S7C.pdf", width=6, height=4,dpi=320) 

######## panel D
print("Making UMAP by sex")
DimPlot(integrated_samples, reduction = "umap", label = FALSE, repel = TRUE, raster = FALSE,
		pt.size=0.1, group.by="sex", shuffle=TRUE, cols = c("M" = sex_pal[1], "F" = sex_pal[2])) + 
	labs(title="") + NoAxes() + NoLegend() 
ggsave(file="../../../paper_figures/S7D.pdf", width=5, height=5,dpi=320)

######## panel E
nuclei_counts <- data.frame(table(integrated_samples@meta.data$celltype_names, 
                                      integrated_samples@meta.data$sampleID, 
                                      integrated_samples@meta.data$sex, 
                                      integrated_samples@meta.data$condition))
colnames(nuclei_counts) <- c("cluster", "sampleID", "sex", "condition", "count") 

a <- ggplot(integrated_samples@meta.data, aes(fill=sex, x=sex)) + 
    ylab("# of nuclei") + scale_fill_manual(values=sex_pal, breaks=c('M', 'F'), labels=c('male', 'female')) +
    geom_bar(show.legend=F) + theme_classic() +
    theme(legend.text = element_text(size=15), legend.title = element_text(size=15),
          axis.title.x = element_text(size=15, color="black"), 
          axis.text.x = element_text(size=15, color="black"),
          axis.title.y = element_text(size=15, color="black"), 
          axis.text.y = element_text(size=15, color="black"))  +
    scale_y_continuous(breaks=c(0,20000,40000,60000), 
                       labels=c("0", "20K","40K","60K")) + 
    scale_x_discrete(limits=c("M","F"))
b <- ggplot(nuclei_counts, aes(fill=sex, y=count, x=cluster)) + 
    ylab("% of nuclei") + scale_fill_manual(values=sex_pal, breaks=c('M', 'F'), labels=c('male', 'female'))  +
    geom_bar(position="fill", stat="identity", show.legend=F) +  coord_flip() + 
    theme_bw() + geom_hline(yintercept=0.5, linetype="dashed") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background=element_rect(colour="black"),
          legend.text = element_text(size=15), legend.title = element_text(size=15),
          axis.title.x = element_text(size=15, vjust=-0.25), 
          axis.text.x = element_text(size=15, color="black"),
          axis.title.y = element_blank(), 
          axis.text.y = element_text(size=15, color="black")) + 
    scale_x_discrete(limits=rev(c("neuron","astro","endo","micro","OPC","oligo"))) +
    scale_y_continuous(labels = function(x) format(x*100, digits=2, nsmall=0))
qc_layout <- c(patchwork::area(t=1,b=4,l=1,r=3),
               patchwork::area(t=1,b=4,l=4,r=10))
plot_all <- a + b + patchwork::plot_layout(design=qc_layout)
ggsave("../../../paper_figures/S7E.pdf", width=6, height=2,dpi=320)    

######## panel F
a <- ggplot(data=integrated_samples@meta.data, aes(x=sex, y=nCount_RNA, fill=sex)) + 
    geom_violin(show.legend=F, trim=F, lwd=0.8) + NoLegend() +
    scale_fill_manual(values=sex_pal, breaks=c('M', 'F'), labels=c('male', 'female')) +
    geom_boxplot(width=0.05, fill="white", outlier.size=0, lwd=0.8) + 
    theme_classic() + ylab("# UMIs") + 
    theme(axis.text.x = element_text(size=15, color="black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=15, color="black"),
        axis.title.y = element_text(size=15, color="black")) +
    scale_x_discrete(limits=c("M","F")) 
b <- ggplot(data=integrated_samples@meta.data, aes(x=sex, y=nFeature_RNA, fill=sex)) + 
    geom_violin(show.legend=F, trim=F, lwd=0.8) + NoLegend() +
    scale_fill_manual(values=sex_pal, breaks=c('M', 'F'), labels=c('male', 'female')) +
    geom_boxplot(width=0.05, fill="white", outlier.size=0, lwd=0.8) + 
    theme_classic() + ylab("# genes") + 
    theme(axis.text.x = element_text(size=15, color="black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=15, color="black"),
        axis.title.y = element_text(size=15, color="black")) +
    scale_x_discrete(limits=c("M","F"))
a | b 
ggsave("../../../paper_figures/S7F.pdf", width=6, height=2,dpi=320)    

######## panel G
print("Making UMAP by sex")
DimPlot(integrated_samples, reduction = "umap", label = FALSE, repel = TRUE, raster = FALSE,
		pt.size=0.1, group.by="condition", shuffle=TRUE, 
		cols = c("control"="gray", "sated"=hungry_sated_pal[1], "hungry"=hungry_sated_pal[2])) + 
	labs(title="") + NoAxes() + NoLegend() 
ggsave(file="../../../paper_figures/S7G.pdf", width=5, height=5,dpi=320)

######## panel H
a <- ggplot(data=integrated_samples@meta.data, aes(x=condition, y=nCount_RNA, fill=condition)) + 
    geom_violin(show.legend=FALSE, trim=F, lwd=0.8) + NoLegend() +
    scale_fill_manual(values=c("gray",hungry_sated_pal), breaks=c('control', 'sated', "hungry")) +
    geom_boxplot(width=0.05, fill="white", outlier.size=0, lwd=0.8) + 
    theme_classic() + ylab("# UMIs") + 
    theme(axis.text.x = element_text(size=15, color="black", angle=30, hjust=1, vjust=1),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=15, color="black"),
          axis.title.y = element_text(size=15, color="black")) +
    scale_y_continuous(breaks=c(0,2000,4000,6000), labels=c("0", "2K","4K","6K")) + 
    scale_x_discrete(limits=c('control', 'sated', "hungry")) 
b <- ggplot(data=integrated_samples@meta.data, aes(x=condition, y=nFeature_RNA, fill=condition)) + 
    geom_violin(show.legend=FALSE, trim=F, lwd=0.8) + NoLegend() +
    scale_fill_manual(values=c("gray",hungry_sated_pal), breaks=c('control', 'sated', "hungry")) +
    geom_boxplot(width=0.05, fill="white", outlier.size=0, lwd=0.8) + 
    theme_classic() + ylab("# genes") + 
    theme(axis.text.x = element_text(size=15, color="black", angle=30, hjust=1, vjust=1),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=15, color="black"),
          axis.title.y = element_text(size=15, color="black")) +
    scale_y_continuous(breaks=c(0,1000,2000), labels=c("0", "1K","2K")) + 
    scale_x_discrete(limits=c('control', 'sated', "hungry"))
a | b 
ggsave("../../../paper_figures/S7H.pdf", width=6, height=2.5,dpi=320)    

#####################

print("Reading in Seurat objects split by celltype") 
all_da <- readRDS(file="../../../seurat_objects/all_da.RDS")  
all_gaba <- readRDS(file="../../../seurat_objects/all_gaba.RDS") 
all_glut <- readRDS(file="../../../seurat_objects/all_glut.RDS") 
celltypes <- c("all_da", "all_gaba", "all_glut") 

# List of genes relevant to feeding or social hormone receptors
food_receptors <- c("Cckar", "Cckbr","Ghsr", "Gipr", "Gpr171", "Hcrtr1", "Hcrtr2", "Insr",  "Lepr", 
                    "Mc1r", "Mc3r", "Mc4r", "Npy1r", "Pparg", "Slc16a2", "Thra", "Thrb", "Trhr", "Tshr") 
food_enzymes <- c("Pcsk1", "Pcsk2", "Cpe")
social_receptors <- c("Ar", "Avpr1a", "Esr1", "Esr2",  "Kiss1r", "Oxtr", "Pgr", "Prlr")  
social_enzymes <- c("Cyp11a1", "Hsd3b2", "Srd5a1", "Srd5a3")
food_genes <- c(food_receptors, food_enzymes)
social_genes <- c(social_receptors, social_enzymes)
food_social_genes <- c(food_genes, social_genes)

######## panel I
print("Making barplot showing percentage of nuclei expressing each food/social gene in each neuron subtype")
expr_df <- read.csv(file="../../../analysis/allConditions_expression_df_across_celltypes.csv")
categories <- c(rep("feeding receptor",length(food_receptors)), rep("feeding enzyme",length(food_enzymes)), 
                rep("social receptor",length(social_receptors)), rep("social enzyme",length(social_enzymes))) 
expr_df$category <- factor(categories, levels=c("feeding receptor", "feeding enzyme",
                                                "social receptor", "social enzyme"))

cols <- c("all_da"="#65A286", "all_gaba"="#965cb5", "all_glut"="#e8d568")
ggplot(expr_df, aes(gene, (pct_cells*100), fill=celltype)) +
    geom_col(position="dodge") +
    theme_classic() + 
    facet_grid(cols=vars(category), scales="free_x", space="free_x", labeller = label_wrap_gen(width=10)) +
    theme(strip.text = element_text(color="black", size=15),
          strip.background = element_rect(fill="#EBEBEB", color=NA), 
          panel.spacing= unit(1,"lines"),
          plot.title = element_text(size=16, color="black", hjust=0.5, vjust=0.5),
          plot.margin = margin(10,100,10,10,unit="pt"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=16), 
          axis.text.x = element_text(color="black", size=12, angle=45, hjust=1, vjust=1),
          axis.title.y = element_text(size=16), 
          axis.text.y = element_text(size=12, color="black"),
          legend.position = c(1.05,0.6), 
          legend.spacing.x = unit(0.1, 'cm'),
          legend.spacing.y = unit(0.3, 'cm'),
          legend.title = element_blank(),
          legend.text = element_text(size=15, color="black", margin=margin(r=30,unit="pt"))) +
    scale_y_continuous(limits=c(0,75), expand=c(0.01,0)) + 
    guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=cols, labels=c("dopamine","gaba","glut")) +
    labs(y="% nuclei", title="expression level across cell types") 
ggsave(file="../../../paper_figures/S7I.pdf", width=15, height=5,dpi=320)

######## panel J
print("Making heatmap showing significance of overlap for each food/social pair in gaba nuclei")
g_pval_mat <- readRDS(file="../../../analysis/g_pval_mat.RDS")
gaba_gpval_heatmap <- pheatmap(g_pval_mat[,,2], main="gaba", cluster_cols=FALSE, cluster_rows=FALSE, 
                             fontsize=12, angle_col=c("45"), border_color="#ededed", breaks=c(0,0.001,0.01,0.05,1), cellwidth=18,
                             color=colorRampPalette(rev(brewer.pal(n=9,name="YlOrRd")))(5), gaps_row=c(8), gaps_col=c(19))
save_pheatmap_pdf(gaba_gpval_heatmap, 
                  "../../../paper_figures/S7J.pdf", 7, 4.5)

######## panel K
g_pval_mat <- readRDS(file="../../../analysis/g_pval_mat.RDS")
print("Making heatmap showing significance of overlap for each food/social pair in glutamate nuclei")
glut_gpval_heatmap <- pheatmap(g_pval_mat[,,3], main="glut", cluster_cols=FALSE, cluster_rows=FALSE, 
                             fontsize=12, angle_col=c("45"), border_color="#ededed", breaks=c(0,0.001,0.01,0.05,1), cellwidth=18,
                             color=colorRampPalette(rev(brewer.pal(n=9,name="YlOrRd")))(5), gaps_row=c(8), gaps_col=c(19))
save_pheatmap_pdf(glut_gpval_heatmap, 
                  "../../../paper_figures/S7K.pdf", 7, 4.5)

print("All done!")





