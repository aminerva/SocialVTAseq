# Script for creating the plots in Figure 4

print("Loading packages")
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("Matrix"))
suppressPackageStartupMessages(library("NatParksPalettes"))         
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("tidyr"))
source("./save_heatmap.R")
set.seed(305)

pcs <- 10 
res <- 0.3 

# Define color palettes
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
food_social_pal <- natparks.pals("KingsCanyon", 6)
food_social_pal <- c(food_social_pal[2], # food
                     food_social_pal[3], # food and social
                     food_social_pal[4], # social
                     food_social_pal[6]) # none

print("Reading in csvs with previously calculated DEGs")
da_degs <- read.csv(file="../../../analysis/da_degs.csv")
gaba_degs <- read.csv(file="../../../analysis/gaba_degs.csv")
glut_degs <- read.csv(file="../../../analysis/glut_degs.csv")
all_hungry_sated_degs <- read.csv(file="../../../analysis/all_hungry_sated_degs.csv")

######## panel B
print("Making volcano plot for each celltype")
cols <- c("down"=hungry_sated_pal[3],"no"=hungry_sated_pal[2],"up"=hungry_sated_pal[1])
options(repr.plot.width = 15, repr.plot.height = 5) 
a <- ggplot(da_degs, aes(x=avg_log2FC, y=-log10(p_val_adj), color=diffexpressed)) +
        geom_point(size = 2) + xlim(-1,1.2) + 
        scale_color_manual(values=cols) + 
        ggtitle("dopamine", "adjusted p<0.01 and log2FC>0.38") +
        labs(x="higher in hungry <--   log2FC   --> higher in sated", y="-10log(adj.p)") +
        NoLegend() + scale_y_continuous(expand=c(0,0), limits=c(-5,200)) +
        theme(axis.text.x = element_text(size=10, color="black"),
              axis.text.y = element_text(size=10, color="black"),
              axis.title.x = element_text(size=12, color="black"),
              axis.title.y = element_text(size=12, color="black"),
              plot.title = element_text(size=12))  
b <- ggplot(gaba_degs, aes(x=avg_log2FC, y=-log10(p_val_adj), color=diffexpressed)) +
        geom_point(size = 2) + xlim(-1,1.2) + 
        scale_color_manual(values=cols) +
        ggtitle("gaba", "adjusted p<0.01 and log2FC>0.38") +
        labs(x="higher in hungry <--   log2FC   --> higher in sated", y="-10log(adj.p)") +
        NoLegend() + scale_y_continuous(expand=c(0,0), limits=c(-5,200)) +
        theme(axis.text.x = element_text(size=10, color="black"),
              axis.text.y = element_text(size=10, color="black"),
              axis.title.x = element_text(size=12, color="black"),
              axis.title.y = element_text(size=12, color="black"),
              plot.title = element_text(size=12))  
c <- ggplot(glut_degs, aes(x=avg_log2FC, y=-log10(p_val_adj), color=diffexpressed)) +
        geom_point(size = 2) + xlim(-1,1.2) + 
        scale_color_manual(values=cols) + 
        ggtitle("glut", "adjusted p<0.01 and log2FC>0.38") +
        labs(x="higher in hungry <--   log2FC   --> higher in sated", y="-10log(adj.p)") +
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
ggsave(file="../../../paper_figures/4B_extra_barplot.pdf",
       width=3.5, height=4,dpi=320)

######## panels C and D made using imaging data and python code 

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
social_enzymes <- c("Cyp11a1", "Hsd3b2", "Srd5a1", "Srd5a3") # "Cyp19a1", "Srd5a2", 
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
                         midpoint = 0, guide = "colourbar", aesthetics="fill") +
    theme(axis.text.y = element_blank()) 
ggsave("../../../paper_figures/4F.pdf", width=10, height=3,dpi=320)

######## panel G
gene1 <- "Insr"
print(paste("Making feature plot for", gene1, "in dopamine neurons"))
if (gene1 %in% food_genes) {
        fp <- FeaturePlot(all_da, features=gene1, slot="data", 
                     order=TRUE, cols=c("#e8e8e8","#855e41"), 
                         pt.size=0.5) + labs(title=gene1) + NoAxes()
    } else if (food_social_genes[gene] %in% social_genes) {
        fp <- FeaturePlot(all_da, features=gene1, slot="data", 
                     order=TRUE, cols=c("#e8e8e8","#8BB1DD"), 
                         pt.size=0.5) + labs(title=gene1) + NoAxes()
    } 
savename1 <- paste0("../../../paper_figures/4G.pdf")
    ggsave(savename1, width=4, height=3,dpi=320)

######## panel H
gene2 <- "Ar"
print(paste("Making feature plot for", gene2, "in dopamine neurons"))
if (gene2 %in% food_genes) {
        fp <- FeaturePlot(all_da, features=gene2, slot="data", 
                     order=TRUE, cols=c("#e8e8e8","#855e41"), 
                         pt.size=0.5) + labs(title=gene2) + NoAxes()
} else if (gene2 %in% social_genes) {
    fp <- FeaturePlot(all_da, features=gene2, slot="data", 
                 order=TRUE, cols=c("#e8e8e8","#8BB1DD"), 
                     pt.size=0.5) + labs(title=gene2) + NoAxes()
} 
savename2 <- paste0("../../../paper_figures/4H.pdf")
    ggsave(savename2, width=4, height=3,dpi=320)

######## panel I
print(paste("Making overlap feature plot for", gene1, "and", gene2, "in dopamine neurons"))
cells <- WhichCells(all_da, expression = Insr > 1 & Ar > 1)
DimPlot(all_da, cells.highlight=cells, cols="#e8e8e8", cols.highlight="#F2C27B",  
        pt.size=0.8, sizes.highlight=0.8) + 
    NoAxes() + NoLegend()
ggsave(file="../../../paper_figures/4I.pdf", width=5, height=3, dpi=320)

######## panel J
gene1 <- "Insr"
gene2 <- "Ar"
print(paste("Making pie chart for", gene1, "and", gene2, "in dopamine neurons"))
pie_ct <- readRDS(file="../../../analysis/pie_counts.RDS")$all_da[[85]]
pie_lab <- readRDS(file="../../../analysis/pie_labs.RDS")$all_da[[85]]
pie_savename <- paste0("../../../paper_figures/4J_top.pdf")
pdf(file=pie_savename)
par(lwd=3)
pie(pie_ct, pie_lab, border="white", init.angle=90, clockwise=TRUE, main="dopamine",
    col = alpha(c(food_social_pal[1], food_social_pal[2], food_social_pal[3], food_social_pal[4]), 0.8)) 
dev.off()

print("Making histogram showing simulated null distribution of overlap between Insr and Ar vs real overlap")
ntotal_da <- length(Cells(all_da))
noverlap <- readRDS(file="../../../analysis/real_overlaps.RDS")$all_da[[85]]
simulated <- readRDS(file="../../../analysis/sim_lists.RDS")$all_da[[85]]  
sim_df <- data.frame(Simulated=simulated) %>% 
              pivot_longer(everything(), names_to="method", values_to="overlap") %>%
                  as.data.frame()
sim_pval <- readRDS(file="../../../analysis/greater_simulation_pvals.RDS")$all_da[[85]]
bw=0.05
ggplot(sim_df, aes(x=(overlap/ntotal_da)*100, fill=method)) +
      geom_histogram(binwidth=bw, alpha=0.25, col="white") +
      geom_density(alpha=0, color="gray", aes(y = bw * after_stat(count)), bw=bw, size=1) +
      geom_vline(xintercept=(noverlap/ntotal_da)*100, lty=2, color="#f2c75b", linewidth=1) +
      theme_bw() +
      theme(text=element_text(size=12), legend.position="none", 
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_blank(), panel.background = element_blank(), 
              axis.line = element_line(colour = "black"), 
              axis.title.x = element_text(size=15), 
              axis.text.x = element_text(size=15, color="black"),
              axis.title.y = element_text(size=15), 
              axis.text.y = element_text(size=15, color="black")) +
      scale_fill_manual(values=c("gray")) + labs(x="da nuclei expressing \n both genes (%)", y="count") +
      scale_x_continuous(limits=c(2.2,3.1), breaks=c(2.25,2.5,2.75,3.0), expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) + 
      annotate(geom="text", label=paste0("p=",round(sim_pval,3)), x=2.3, y=1800) + 
      annotate(geom="text", label="true \noverlap", x=2.75, y=1750, color="#f2c75b", hjust=0, lineheight=0.8)
ggsave("../../../paper_figures/4J_bottom.pdf", width=4, height=3,dpi=320) 

######## panel K
g_pval_mat <- readRDS(file="../../../analysis/g_pval_mat.RDS")
da_gpval_heatmap <- pheatmap(g_pval_mat[,,1], main="dopamine", cluster_cols=FALSE, cluster_rows=FALSE, 
                             fontsize=12, angle_col=c("45"), border_color="#ededed", 
                             breaks=c(0,0.001,0.01,0.05,1), cellwidth=18,
                             color=colorRampPalette(rev(brewer.pal(n=9,name="YlOrRd")))(5), 
                             gaps_row=c(8), gaps_col=c(19))
save_pheatmap_pdf(da_gpval_heatmap, 
                  "../../../paper_figures/4K.pdf", 7, 4.5)



######## panel L
print("Making bar plot with percentage of gene pairs overlapping more than chance in each neuron subtype")
g_pval_mat <- readRDS(file="../../../analysis/g_pval_mat.RDS")
all_gpval_df <- readRDS(file="../../../analysis/all_gpval_df.RDS")

cols <- c("all_da"="#65A286", "all_gaba"="#965cb5", "all_glut"="#e8d568")
all_gpval_df[all_gpval_df$celltype %in% c("all_da","all_gaba","all_glut"),] %>% 
    group_by(celltype) %>%
        summarise(n=length(which(pvalue<0.05)), pct=round(n/(nrow(g_pval_mat)*ncol(g_pval_mat))*100, 0)) %>%
            ggplot(mapping=aes(x=celltype,y=pct,fill=celltype)) + 
                geom_col() + scale_fill_manual(values=cols) +
                theme_bw() + NoLegend() + 
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(), 
                      panel.border = element_blank(),
                      plot.title = element_text(hjust=0.5, size=15, color="black"),
                      axis.line = element_line(colour="black"),
                      axis.title.x = element_text(size=15), 
                      axis.text.x = element_text(size=15, color="black", angle=0),
                      axis.title.y = element_text(size=18, color="black"), 
                      axis.text.y = element_text(size=15, color="black")) +
                scale_y_continuous(expand=expansion(mult=c(0,0.05)), breaks=c(0,5,10,15,20)) +
                labs(y="percentage of pairs \n overlapping more than chance", x=" ") +
                scale_x_discrete(expand=c(0.2,0.2), labels=c("da","gaba","glut")) +
                annotate(geom="text", label="*", x=2, y=19.4, size=8) +
                annotate(geom="segment", y=19.3, yend=19.3, x=1, xend=3) +
                annotate(geom="text", label="ns", x=1.5, y=21.1, size=5) +
                annotate(geom="segment", y=20.4, yend=20.4, x=1, xend=2)
ggsave(file="../../../paper_figures/4L.pdf",width=3,height=4,dpi=320)

print("All done!")














