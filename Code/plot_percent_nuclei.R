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
suppressPackageStartupMessages(library("sctransform"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("tidyr"))

source("./std_error.R")
source("./snRNA_functions.R")
source("./save_heatmap.R")

set.seed(305)

# Define color palettes 
fs_pal <- natparks.pals("KingsCanyon", 6)
fs_pal <- c(fs_pal[2], # food
            fs_pal[3], # food and social
            fs_pal[4], # social
            fs_pal[6]) # none
hs_pal <- natparks.pals("Glacier", 5)
hs_pal <- c(hs_pal[2], # sated
            hs_pal[4]) # hungry

print("Reading in the Seurat objects that have been split by celltype and condition")
ctrl_da <- readRDS(file="../../../seurat_objects/ctrl_da_obj.RDS") 
ctrl_gaba <- readRDS(file="../../../seurat_objects/ctrl_gaba_obj.RDS") 
ctrl_glut <- readRDS(file="../../../seurat_objects/ctrl_glut_obj.RDS") 

hungry_da <- readRDS(file="../../../seurat_objects/hungry_da_obj.RDS") 
hungry_gaba <- readRDS(file="../../../seurat_objects/hungry_gaba_obj.RDS") 
hungry_glut <- readRDS(file="../../../seurat_objects/hungry_glut_obj.RDS") 

sated_da <- readRDS(file="../../../seurat_objects/sated_da_obj.RDS") 
sated_gaba <- readRDS(file="../../../seurat_objects/sated_gaba_obj.RDS") 
sated_glut <- readRDS(file="../../../seurat_objects/sated_glut_obj.RDS") 

# List of genes relevant to feeding or social hormone receptors
food_receptors <- c("Cckar", "Cckbr","Ghsr", "Gipr", "Gpr171", "Hcrtr1", "Hcrtr2", "Insr",  "Lepr", 
                    "Mc1r", "Mc3r", "Mc4r", "Npy1r", "Pparg", "Slc16a2", "Thra", "Thrb", "Trhr", "Tshr") 
food_enzymes <- c("Pcsk1", "Pcsk2", "Cpe")
social_receptors <- c("Ar", "Avpr1a", "Esr1", "Esr2",  "Kiss1r", "Oxtr", "Pgr", "Prlr")  
social_enzymes <- c("Cyp19a1", "Cyp11a1", "Hsd3b2", "Srd5a1", "Srd5a2", "Srd5a3") 
food_genes <- c(food_receptors, food_enzymes)
social_genes <- c(social_receptors, social_enzymes)
food_social_genes <- c(food_genes, social_genes)

typenames <- c("ctrl_da",   "ctrl_gaba",   "ctrl_glut", 
               "hungry_da", "hungry_gaba", "hungry_glut",
               "sated_da",  "sated_gaba",  "sated_glut")

# Make a df containing expression info for all nuclei across different celltypes 
print("Calculating gene expression metrics for dataframe")
expr_df <- calcGeneExpressionMetrics(celltype_names=typenames,
                                     genes_of_interest=c(food_receptors,food_enzymes,
                                                         social_receptors,social_enzymes))
# Add a column indicating category of each gene
categories <- c(rep("feeding receptor",length(food_receptors)), rep("feeding enzyme",length(food_enzymes)), 
                rep("social receptor",length(social_receptors)), rep("social enzyme",length(social_enzymes))) 
expr_df$category <- rep(categories, length(typenames))
expr_df$category <- factor(categories, levels=c("feeding receptor", "feeding enzyme",
                                           "social receptor", "social enzyme"))
write.csv(expr_df, file="../../../analysis/all_expression_df.csv")

#######################

print("Plotting baseline expression across celltypes")
options(repr.plot.width = 15, repr.plot.height = 5) 
filt <- filter(expr_df, 
               celltype %in% c("ctrl_da","ctrl_gaba","ctrl_glut"))[c("gene","celltype","pct_cells","category")] 

p <- ggplot(filt, aes(gene, (pct_cells*100), fill=celltype)) +
        geom_col(position="dodge") +
        theme_classic() + 
        facet_grid(cols=vars(category), scales="free_x", space="free_x", labeller=label_wrap_gen(width=10)) +
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
        guides(fill=guide_legend(ncol=1)) + scale_fill_discrete(labels=c("dopamine","gaba","glut")) +
        labs(y="% nuclei", title="expression level across cell types") 
p
ggsave(file="../../../plots/control_food_social_expression.pdf",
       width=15, height=5,dpi=320)
ggsave(file="../../../plots/control_food_social_expression.png",
       width=15, height=5,dpi=320)

#######################

print("Plotting expression in DA nuclei across hunger states")
options(repr.plot.width = 15, repr.plot.height = 5) 
filt <- filter(expr_df, celltype %in% c("ctrl_da",
                                        "hungry_da",
                                        "sated_da"))[c("gene","celltype","pct_cells","category")] 
p <- ggplot(filt, aes(gene, (pct_cells*100), fill=celltype)) +
        geom_col(position="dodge") +
        theme_classic() + 
        facet_grid(cols=vars(category), scales="free_x", space="free_x", labeller=label_wrap_gen(width=10)) +
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
        guides(fill=guide_legend(ncol=1)) + 
        scale_fill_manual(labels=c("control da","hungry da","sated da"), 
                          values=c("gray", hs_pal[2], hs_pal[1])) +
        labs(y="% nuclei", title="expression level states (da nuclei)") 
p
ggsave(file="../../../plots/da_food_social_expression_across_states.pdf",
       width=15, height=5,dpi=320)
ggsave(file="../../../plots/da_food_social_expression_across_states.png",
       width=15, height=5,dpi=320)

#######################

print("Plotting expression in GABA nuclei across hunger states")
options(repr.plot.width = 15, repr.plot.height = 5) 
filt <- filter(expr_df, celltype %in% c("ctrl_gaba",
                                        "hungry_gaba",
                                        "sated_gaba"))[c("gene","celltype","pct_cells","category")] 
p <- ggplot(filt, aes(gene, (pct_cells*100), fill=celltype)) +
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
        guides(fill=guide_legend(ncol=1)) + 
        scale_fill_manual(labels=c("control gaba","hungry gaba","sated gaba"), 
                          values=c("gray", hs_pal[2], hs_pal[1])) +
        labs(y="% nuclei", title="expression level across states (gaba nuclei)") 
p
ggsave(file="../../../plots/gaba_food_social_expression_across_states.pdf",
       width=15, height=5,dpi=320)
ggsave(file="../../../plots/gaba_food_social_expression_across_states.png",
       width=15, height=5,dpi=320)

#######################

print("Plotting expression in glut nuclei across hunger states")
options(repr.plot.width = 15, repr.plot.height = 5) 
filt <- filter(expr_df, celltype %in% c("ctrl_glut",
                                        "hungry_glut",
                                        "sated_glut"))[c("gene","celltype","pct_cells","category")] 
p <- ggplot(filt, aes(gene, (pct_cells*100), fill=celltype)) +
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
        guides(fill=guide_legend(ncol=1)) + 
        scale_fill_manual(labels=c("control glut","hungry glut","sated glut"), 
                          values=c("gray", hs_pal[2], hs_pal[1])) +
        labs(y="% nuclei", title="expression level across states (glut nuclei)") 
p
ggsave(file="../../../plots/glut_food_social_expression_across_states.pdf",
       width=15, height=5,dpi=320)
ggsave(file="../../../plots/glut_food_social_expression_across_states.png",
       width=15, height=5,dpi=320)

print("All done!")







