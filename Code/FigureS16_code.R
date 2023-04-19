# Script for creating the plots in Figure S16

print("Loading packages")
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2")) 
suppressPackageStartupMessages(library("NatParksPalettes"))         
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("RColorBrewer"))    
suppressPackageStartupMessages(library("Seurat"))
source("./save_heatmap.R")
set.seed(305)

print("Reading in Seurat objects") 
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

######## panel A
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
ggsave(file="../../../paper_figures/S16A.pdf", width=15, height=5,dpi=320)

######## panel B
print("Making heatmap showing significance of overlap for each food/social pair in gaba nuclei")
g_pval_mat <- readRDS(file="../../../analysis/g_pval_mat.RDS")
gaba_gpval_heatmap <- pheatmap(g_pval_mat[,,2], main="gaba", cluster_cols=FALSE, cluster_rows=FALSE, 
                             fontsize=12, angle_col=c("45"), border_color="#ededed", breaks=c(0,0.001,0.01,0.05,1), cellwidth=18,
                             color=colorRampPalette(rev(brewer.pal(n=9,name="YlOrRd")))(5), gaps_row=c(8), gaps_col=c(19))
save_pheatmap_pdf(gaba_gpval_heatmap, 
                  "../../../paper_figures/S16B.pdf", 7, 4.5)

######## panel C
g_pval_mat <- readRDS(file="../../../analysis/g_pval_mat.RDS")
print("Making heatmap showing significance of overlap for each food/social pair in glutamate nuclei")
glut_gpval_heatmap <- pheatmap(g_pval_mat[,,3], main="glut", cluster_cols=FALSE, cluster_rows=FALSE, 
                             fontsize=12, angle_col=c("45"), border_color="#ededed", breaks=c(0,0.001,0.01,0.05,1), cellwidth=18,
                             color=colorRampPalette(rev(brewer.pal(n=9,name="YlOrRd")))(5), gaps_row=c(8), gaps_col=c(19))
save_pheatmap_pdf(glut_gpval_heatmap, 
                  "../../../paper_figures/S16C.pdf", 7, 4.5)

print("All done!")






