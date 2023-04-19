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
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("reldist"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("tidyr"))

source("./std_error.R")
source("./snRNA_functions.R")
source("./save_heatmap.R")

set.seed(305)

# Define color palettes 
food_social_pal <- natparks.pals("KingsCanyon", 6)
food_social_pal <- c(food_social_pal[2], # food
                     food_social_pal[3], # food and social
                     food_social_pal[4], # social
                     food_social_pal[6]) # none
hungry_sated_pal <- natparks.pals("Glacier", 5)
hungry_sated_pal <- c(hungry_sated_pal[2], # sated
                      hungry_sated_pal[4]) # hungry

print("Reading in seurat objects")
all_da <- readRDS(file="../../../seurat_objects/all_da.RDS") 
all_gaba <- readRDS(file="../../../seurat_objects/all_gaba.RDS") 
all_glut <- readRDS(file="../../../seurat_objects/all_glut.RDS") 
celltypes <- c("all_da", "all_gaba", "all_glut") 

# List of genes relevant to feeding or social hormone receptors
food_receptors <- c("Cckar", "Cckbr","Ghsr", "Gipr", "Gpr171", "Hcrtr1", "Hcrtr2", "Insr",  "Lepr", 
                    "Mc1r", "Mc3r", "Mc4r", "Npy1r", "Pparg", "Slc16a2", "Thra", "Thrb", "Trhr", "Tshr") 
food_enzymes <- c("Pcsk1", "Pcsk2", "Cpe")
social_receptors <- c("Ar", "Avpr1a", "Esr1", "Esr2",  "Kiss1r", "Oxtr", "Pgr", "Prlr")  
social_enzymes <- c("Cyp11a1", "Hsd3b2", "Srd5a1", "Srd5a3") # "Cyp19a1", "Srd5a2",
food_genes <- c(food_receptors, food_enzymes)
social_genes <- c(social_receptors, social_enzymes)
food_social_genes <- c(food_genes, social_genes)

print("Calculating overlap and plotting")
ID_list <- readRDS(file="../../../analysis/food_social_nuclei_ID_list_allConditions.RDS")

pie_counts <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list())
pie_percents <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list())
pie_labs <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list())
real_overlaps <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list()) # true % of coexpressing cells
sim_lists <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list()) # num cells expr both genes from one simulation
g_sim_pvals <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list()) # pvals for > chance overlap
l_sim_pvals <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list()) # pvals for < chance overlap

for (type in seq_along(celltypes)) { # Iterate through the three cell types
    
    print(paste("Starting calculations for", celltypes[type]))
    
    for (i in seq_along(food_genes)) { # Iterate through food genes
        
        for (j in seq_along(social_genes)) { # Iterate through social genes


            celltype <- celltypes[type] 
            food_gene <- food_genes[i]
            social_gene <- social_genes[j]
            ntotal <- 0 
            nfood_gene <- 0
            nsocial_gene <- 0
            nboth <- 0
            ntotal <- ncol(get(celltypes[type]))

            print(paste("Calculating overlap for", food_gene, "and", social_gene))
            
            # Calculate # of cells of given type expressing the food gene from the pair
            nfood_gene <- nfood_gene + length(ID_list[[celltypes[type]]][[food_gene]])

            # Calculate # of cells of given celltype expressing the social gene from the pair
            nsocial_gene <- nsocial_gene + length(ID_list[[celltypes[type]]][[social_gene]])

            # Calculate # of cells coexpressing both genes by finding the length of the overlap of cell IDs
            nboth <- nboth + length(intersect(ID_list[[celltypes[type]]][[food_gene]],
                                              ID_list[[celltypes[type]]][[social_gene]])) 

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
            
            print("Plotting and saving pie chart")
            options(repr.plot.width = 6, repr.plot.height = 5)
            pie_savename <- paste0("../../../plots/",food_gene,"_",social_gene,"_",celltype,"_pie.pdf")
            pdf(file=pie_savename)
            par(lwd=3)
            pie(pie_ct, pie_lab, border="white", init.angle=90, clockwise=TRUE, main=celltype,
                col = alpha(c(food_social_pal[1], food_social_pal[2], food_social_pal[3], food_social_pal[4]), 0.8)) 
            dev.off()

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

            print("Plotting and saving simulated overlap")
            sim_lists[[type]] <- c(sim_lists[[type]], list(simulated)) # Save simulated co-expr to list
            sim_df <- data.frame(Simulated=simulated) %>% 
                      pivot_longer(everything(), names_to="Method", values_to="Overlap") %>%
                          as.data.frame()

            greater_pval <- (sum(simulated >= nboth) + 1) / (10000 + 1) # <0.05 if real >= sim
            g_sim_pvals[[type]] <- c(g_sim_pvals[[type]], list(greater_pval)) 
            less_pval <- (sum(simulated <= nboth) + 1) / (10000 + 1) # <0.05 if real <= sim
            l_sim_pvals[[type]] <- c(l_sim_pvals[[type]], list(less_pval)) 

            print(paste(food_gene, "-", social_gene, "overlap greater than chance pval=", greater_pval))
            print(paste(food_gene, "-", social_gene, "overlap less than chance pval=", less_pval))

            ggplot(sim_df, aes(x=Overlap, fill=Method)) +
                geom_histogram(binwidth=2, col="white") +
                geom_vline(xintercept=nboth, lty=2, color="#F1EE88", size=2) +
                theme_bw() + 
                theme(text=element_text(size=12), legend.position="none", 
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.border = element_blank(), panel.background = element_blank(), 
                      axis.line = element_line(colour = "black")) +
                scale_fill_manual(values=c("gray")) + ylab("Count") +
                labs(title = paste(food_gene,"-",social_gene), 
                     subtitle = paste0("p>real=", round(greater_pval,5), " / p<real=", round(less_pval,5)))
            savename <- paste0("../../../plots/",food_gene,"_",social_gene,"_",celltype,"_overlap_simulation.pdf")
            ggsave(file=savename,width=6,height=5,dpi=320)
        }

        print(paste("All social gene comparisons with", food_genes[i], "for", celltypes[type], "nuclei done!"))   
    }
    print(paste(celltypes[type], "all done!"))
}



