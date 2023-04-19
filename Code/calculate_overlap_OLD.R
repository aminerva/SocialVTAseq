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
suppressPackageStartupMessages(library("purrr"))
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

print("Calculating metrics for those nuclei that express each food/social gene")
fs_df <- calcExpressingNuclei(celltype_names=typenames,
                                genes_of_interest=c(food_receptors,food_enzymes,
                                                    social_receptors,social_enzymes))
# Add a column indicating category of each gene
categories <- c()
for (gene in fs_df$gene){
    if (gene %in% food_receptors) {
        cat <- c("feeding receptor")
    } else if (gene %in% food_enzymes) {
        cat <- c("feeding enzyme") 
    } else if (gene %in% social_receptors) {
        cat <- c("social receptor")
    } else if (gene %in% social_enzymes) {
        cat <- c("social enzyme")
    }
    categories <- c(categories, cat)
}
fs_df$category <- categories
write.csv(fs_df, file="../../../analysis/food_social_expressing_nuclei.csv")

print("Grabbing the IDs of the nuclei that express each gene")
ID_list <- list("ctrl_da"=list(), "ctrl_gaba"=list(), "ctrl_glut"=list(),
                "hungry_da"=list(), "hungry_gaba"=list(), "hungry_glut"=list(),
                "sated_da"=list(), "sated_gaba"=list(), "sated_glut"=list())
for (type in 1:length(typenames)) { #Iterate through the cell types
    for (i in 1:length(food_social_genes)) { # Iterate through all genes
        # Grab nucleus IDs (these are the nuclei within that celltype that express a given gene)
        IDs <- fs_df[fs_df$celltype==typenames[type] & fs_df$gene==food_social_genes[i],]$nucleusID
        ID_list[[type]][[food_social_genes[i]]] = c(IDs)
    }
}
saveRDS(ID_list, file="../../../analysis/food_social_nuclei_ID_list.RDS")


pie_counts <- list("da"=list(), "gaba"=list(), "glut"=list())
pie_percents <- list("da"=list(), "gaba"=list(), "glut"=list())
pie_labs <- list("dat_or_th"=list(), "gaba"=list(), "glut"=list())
real_overlaps <- list("da"=list(), "gaba"=list(), "glut"=list()) # true % of coexpressing cells
sim_lists <- list("da"=list(), "gaba"=list(), "glut"=list()) # num cells expr both genes from one simulation
g_sim_pvals <- list("da"=list(), "gaba"=list(), "glut"=list()) # pvals for > chance overlap
l_sim_pvals <- list("da"=list(), "gaba"=list(), "glut"=list()) # pvals for < chance overlap

master_celltypes <- c("da","gaba","glut")
states <- c("ctrl","hungry","sated")
for (type in seq_along(master_celltypes)) { # Iterate through the three cell types
    
    print(paste("Starting calculations for", master_celltypes[type]))
    
    for (i in seq_along(food_genes)) { # Iterate through food genes
        
        for (j in seq_along(social_genes)) { # Iterate through social genes

            celltype <- master_celltypes[type] 
            food_gene <- food_genes[i]
            social_gene <- social_genes[j]
            
            # Grab the names of all three seurat objects for each celltype (ctrl, hungry, sated)
            allstates_celltype <- typenames[grepl(master_celltypes[type],typenames)]
            ntotal <- 0 
            nfood_gene <- 0
            nsocial_gene <- 0
            nboth <- 0

            for (state in seq_along(states)) {
                ntotal <- ntotal + ncol(get(allstates_celltype[state]))

                # Calculate # of cells of given type expressing the food gene from the pair
                nfood_gene <- nfood_gene + length(ID_list[[allstates_celltype[state]]][[food_gene]])

                # Calculate # of cells of given celltype expressing the social gene from the pair
                nsocial_gene <- nsocial_gene + length(ID_list[[allstates_celltype[state]]][[social_gene]])

                # Calculate # of cells coexpressing both genes by finding the length of the overlap of cell IDs
                nboth <- nboth + length(intersect(ID_list[[allstates_celltype[state]]][[food_gene]],
                                                  ID_list[[allstates_celltype[state]]][[social_gene]])) 
            }

            # Add real overlaps to list
            real_overlaps[[type]] <- c(real_overlaps[[type]], nboth) 
            
            # Calculate # of cells expressing only one of the two genes
            nfood_gene_only <- nfood_gene - nboth
            nsocial_gene_only <- nsocial_gene - nboth
            
            # Calculate # of cells expressing neither gene in the pair
            nnone <- ntotal - nboth - nfood_gene_only - nsocial_gene_only

            ## PIE CHARTS ##

            # Put numbers into format for plotting in pie chart
            pie_ct <- c(nfood_gene_only, nboth, nsocial_gene_only, nnone) 
            pie_counts[[celltype]] <- c(pie_counts[[celltype]], list(pie_ct)) # Save to list

            # Turn pie counts into percent of total cells of that type
            pie_pct <- round(pie_ct/ntotal*100, 0) 
            pie_percents[[celltype]] <- c(pie_percents[[celltype]], list(pie_pct)) # Save to list

            # Pie labels need to be in same order as data (gene 1 only, overlap, gene 2 only, none)
            pie_lab <- c(paste0(food_gene," (",pie_pct[1],"%)"), 
                         paste0("Both (",pie_pct[2],"%)"), 
                         paste0(social_gene," (",pie_pct[3],"%)"), 
                         paste0("None (",pie_pct[4],"%)"))
            pie_labs[[celltype]] <- c(pie_labs[[celltype]], list(pie_lab)) # Save to list
            

            ## SIMULATED OVERLAP DISTRIBUTION ##

            # Number of nuclei expressing the food gene or the social gene, regardless of coexpression
            sets <- c(nfood_gene, nsocial_gene) 

            # Simulate null distribution of overlap 10000x, given the number of nuclei exprerssing each gene
            simulated <- purrr::map_dbl(seq_len(10000), function(x) {

                # For each # of nuclei expressing the food gene and those expressing the social gene, 
                # Take a sample from ntotal of that size 
                    # The output of this is two vectors, 
                        # which each contain simulated indices/IDs of a cell expressing that gene
                        # The length of each of these vectors is the same as the length of nfood_gene and nsocial_gene
                sim <- map(sets, ~sample(ntotal, .x))

                # Get the intersection of the two simulated lists and find the length of it
                    # aka how many simulated cells express gene1 AND gene2 (simulated co-expression)
                sim <- length(purrr::reduce(sim, intersect))
                return(sim)
            })

            sim_lists[[type]] <- c(sim_lists[[type]], list(simulated)) # Save simulated co-expr to list
            greater_pval <- (sum(simulated >= nboth) + 1) / (10000 + 1) # <0.05 if real >= sim
            g_sim_pvals[[type]] <- c(g_sim_pvals[[type]], list(greater_pval)) 
            less_pval <- (sum(simulated <= nboth) + 1) / (10000 + 1) # <0.05 if real <= sim
            l_sim_pvals[[type]] <- c(l_sim_pvals[[type]], list(less_pval)) 
        } 
        print(paste("All social gene comparisons with", food_genes[i], "for", master_celltypes[type], "nuclei done!"))   
    }
    print(paste(master_celltypes[type], "all done!"))
}

saveRDS(pie_counts, file="../../../analysis/pie_counts.RDS") 
saveRDS(pie_percents, file="../../../analysis/pie_percents.RDS") 
saveRDS(pie_labs, file="../../../analysis/pie_labs.RDS") 
saveRDS(real_overlaps, file="../../../analysis/real_overlaps.RDS") 
saveRDS(sim_lists, file="../../../analysis/sim_lists.RDS")  
saveRDS(g_sim_pvals, file="../../../analysis/greater_simulation_pvals.RDS")
saveRDS(l_sim_pvals, file="../../../analysis/less_simulation_pvals.RDS") 












