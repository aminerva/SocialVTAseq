print("Loading packages")
suppressPackageStartupMessages(library("abind"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("grid"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("gtools"))
suppressPackageStartupMessages(library("likert"))
suppressPackageStartupMessages(library("Matrix"))
suppressPackageStartupMessages(library("metap"))     
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("reldist"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("tidyr"))

source("./std_error.R")
source("./snRNA_functions.R")

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

print("Reading in Seurat object") 
all_da <- readRDS(file="../../../seurat_objects/all_da_integrated_clustered_2000features_10pcs_0.3res.RDS")
DefaultAssay(all_da) <- "RNA"

# List of genes relevant to feeding or social hormone receptors
food_receptors <- c("Cckar", "Cckbr","Ghsr", "Gipr", "Gpr171", "Hcrtr1", "Hcrtr2", "Insr",  "Lepr", 
                    "Mc1r", "Mc3r", "Mc4r", "Npy1r", "Pparg", "Slc16a2", "Thra", "Thrb", "Trhr", "Tshr") 
food_enzymes <- c("Pcsk1", "Pcsk2", "Cpe")
social_receptors <- c("Ar", "Avpr1a", "Esr1", "Esr2",  "Kiss1r", "Oxtr", "Pgr", "Prlr")  
social_enzymes <- c("Cyp11a1", "Hsd3b2", "Srd5a1", "Srd5a3") # "Cyp19a1", "Srd5a2",
food_genes <- c(food_receptors, food_enzymes)
social_genes <- c(social_receptors, social_enzymes)
food_social_genes <- c(food_genes, social_genes)

print("Calculating coexpression of each food/social gene pair in DA, GABA, and glutamate neurons")

print("Reading in the all neurons seurat object")
all_neurons <- readRDS(file="../../../seurat_objects/all_conditions_all_neurons.RDS")
DefaultAssay(all_neurons) <- "RNA"

print("Pulling out just dopamine neurons")
all_da <- subset(all_neurons, subset = Slc6a3 > 0 | Th > 0 ) 
all_da[["nucleusID"]] <- colnames(all_da)
saveRDS(all_da, file="../../../seurat_objects/all_da.RDS")  

print("Pulling out just GABA neurons")
all_gaba <- subset(all_neurons, subset = Slc32a1>0 & Slc6a3==0 & Th==0 & Slc17a6==0& Slc17a7==0 & Slc17a8==0 | Gad1>0 & Slc6a3==0 & Th==0 & Slc17a6==0 & Slc17a7==0 & Slc17a8==0 | Gad2>0 & Slc6a3==0 & Th==0 & Slc17a6==0 & Slc17a7==0 & Slc17a8==0)
all_gaba[["nucleusID"]] <- colnames(all_gaba)
saveRDS(all_gaba, file="../../../seurat_objects/all_gaba.RDS") 

print("Pulling out just glut neurons")
all_glut <- subset(all_neurons, subset = Slc17a6>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0 | Slc17a7>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0 | Slc17a8>0 & Slc6a3==0 & Th==0 & Slc32a1==0 & Gad1==0 & Gad2==0) 
all_glut[["nucleusID"]] <- colnames(all_glut)
saveRDS(all_glut, file="../../../seurat_objects/all_glut.RDS") 

celltypes <- c("all_da", "all_gaba", "all_glut") 

print("Calculating metrics for those nuclei that express each food/social gene")
fs_df <- calcExpressingNuclei(celltype_names=celltypes, genes_of_interest=c(food_receptors,food_enzymes,
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
write.csv(fs_df, file="../../../analysis/food_social_expressing_nuclei_allConditions.csv")

print("Grabbing the IDs of the nuclei that express each gene")
ID_list <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list())
for (type in 1:length(celltypes)) { 
    for (i in 1:length(food_social_genes)) { 
        # Grab nucleus IDs (these are the nuclei within that celltype that express a given gene)
        IDs <- fs_df[fs_df$celltype==celltypes[type] & fs_df$gene==food_social_genes[i],]$nucleusID
        ID_list[[type]][[food_social_genes[i]]] = c(IDs)
    }
}
saveRDS(ID_list, file="../../../analysis/food_social_nuclei_ID_list_allConditions.RDS")

pie_counts <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list())
pie_percents <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list())
pie_labs <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list())
real_overlaps <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list()) # true % of coexpressing cells
sim_lists <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list()) # num cells expr both genes from one simulation
g_sim_pvals <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list()) # pvals for > chance overlap
l_sim_pvals <- list("all_da"=list(), "all_gaba"=list(), "all_glut"=list()) # pvals for < chance overlap

for (type in seq_along(celltypes)) { 
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

            print("Calculating and saving simulated overlap")
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
        }

        print(paste("All social gene comparisons with", food_genes[i], "for", celltypes[type], "nuclei done!"))   
    }
    print(paste(celltypes[type], "all done!"))
}

print("Saving RDS files with overlap calculations")
saveRDS(pie_counts, file="../../../analysis/pie_counts.RDS") 
saveRDS(pie_percents, file="../../../analysis/pie_percents.RDS") 
saveRDS(pie_labs, file="../../../analysis/pie_labs.RDS") 
saveRDS(real_overlaps, file="../../../analysis/real_overlaps.RDS") 
saveRDS(sim_lists, file="../../../analysis/sim_lists.RDS")  
saveRDS(g_sim_pvals, file="../../../analysis/greater_simulation_pvals.RDS")
saveRDS(l_sim_pvals, file="../../../analysis/less_simulation_pvals.RDS") 

ntotal <- c(length(Cells(all_da)), length(Cells(all_gaba)), length(Cells(all_glut)))

print("Creating matrices with overlap info")
for (type in seq_along(celltypes)) {
    
    ## real percent overlap matrices ##
    coexpr_mat <- matrix(unlist(real_overlaps[[celltypes[[type]]]]), 
                                nrow=length(social_genes), 
                                ncol=length(food_genes)) 
    coexpr_mat <- coexpr_mat / ntotal[type] * 100
    rownames(coexpr_mat) <- social_genes
    colnames(coexpr_mat) <- food_genes
    assign(paste0(celltypes[type],"_coexpr"), coexpr_mat)

    ## greater than pval matrices ##
    greater_pval_mat <- matrix(unlist(g_sim_pvals[[celltypes[[type]]]]), 
                               nrow=length(social_genes), 
                               ncol=length(food_genes)) # Pvals for whether overlap is sig > than chance
    rownames(greater_pval_mat) <- social_genes
    colnames(greater_pval_mat) <- food_genes
    assign(paste0(celltypes[type],"_greater_pvals"), greater_pval_mat)  

    ## less than pval matrices ##
    less_pval_mat <- matrix(unlist(l_sim_pvals[[celltypes[[type]]]]), 
                            nrow=length(social_genes), 
                            ncol=length(food_genes)) # Pvals for whether overlap is sig < than chance
    rownames(less_pval_mat) <- social_genes
    colnames(less_pval_mat) <- food_genes
    assign(paste0(celltypes[type],"_less_pvals"), less_pval_mat)
}
coexpr_mat <- abind(all_da_coexpr, all_gaba_coexpr, all_glut_coexpr, along=3)   
g_pval_mat <- abind(all_da_greater_pvals, all_gaba_greater_pvals, all_glut_greater_pvals, along=3)   
l_pval_mat <- abind(all_da_less_pvals, all_gaba_less_pvals, all_glut_less_pvals, along=3)
saveRDS(coexpr_mat, file="../../../analysis/coexpr_mat.RDS")
saveRDS(g_pval_mat, file="../../../analysis/g_pval_mat.RDS")
saveRDS(l_pval_mat, file="../../../analysis/l_pval_mat.RDS")

print("Reshaping the overlap matrices into dataframes for plotting")
for (type in seq_along(celltypes)) {
    
    coexpr_df <- reshape2::melt(coexpr_mat[,,type])
    colnames(coexpr_df) <- c("gene1", "gene2", "pct_both")
    coexpr_df$celltype <- rep(celltypes[type], nrow(coexpr_df))
    assign(paste0(celltypes[type],"_coexpr_df"), coexpr_df)
    
    g_pval_df <- reshape2::melt(g_pval_mat[,,type])
    colnames(g_pval_df) <- c("gene1", "gene2", "pvalue")
    g_pval_df$celltype <- rep(celltypes[type], nrow(g_pval_df))
    assign(paste0(celltypes[type],"_gpval_df"), g_pval_df)
    
    l_pval_df <- reshape2::melt(l_pval_mat[,,type])
    colnames(l_pval_df) <- c("gene1", "gene2", "pvalue")
    l_pval_df$celltype <- rep(celltypes[type], nrow(l_pval_df))
    assign(paste0(celltypes[type],"_lpval_df"), l_pval_df)
}

all_coexpr_df <- rbind(all_da_coexpr_df[c("pct_both", "celltype")], 
                       all_gaba_coexpr_df[c("pct_both", "celltype")], 
                       all_glut_coexpr_df[c("pct_both", "celltype")])
all_gpval_df <- rbind(all_da_gpval_df[c("pvalue", "celltype")], 
                      all_gaba_gpval_df[c("pvalue", "celltype")], 
                      all_glut_gpval_df[c("pvalue", "celltype")])
all_lpval_df <- rbind(all_da_lpval_df[c("pvalue", "celltype")], 
                      all_gaba_lpval_df[c("pvalue", "celltype")], 
                      all_glut_lpval_df[c("pvalue", "celltype")])
saveRDS(all_coexpr_df, file="../../../analysis/all_coexpr_df.RDS")
saveRDS(all_gpval_df, file="../../../analysis/all_gpval_df.RDS")
saveRDS(all_lpval_df, file="../../../analysis/all_lpval_df.RDS")

print("Performing tests of proportions")
pct_sig <- all_gpval_df[all_gpval_df$celltype %in% c("all_da","all_gaba","all_glut"),] %>% 
                   group_by(celltype) %>% summarise(n=length(which(pvalue<0.05)), 
                                             pct=round(n/(nrow(g_pval_mat)*ncol(g_pval_mat))*100, 2))
n_da <- pct_sig[pct_sig$celltype=="all_da",]$n
n_gaba <- pct_sig[pct_sig$celltype=="all_gaba",]$n
n_glut <- pct_sig[pct_sig$celltype=="all_glut",]$n
total_pairs <- length(food_genes) * length(social_genes)

p_da_gaba <- round(prop.test(x = c(n_da, n_gaba), n = c(total_pairs, total_pairs), correct=TRUE)$p.value, 3)
chi_da_gaba <- round(prop.test(x = c(n_da, n_gaba), n = c(total_pairs, total_pairs), correct=TRUE)$statistic, 3)
print(paste0("DA vs GABA: p=", p_da_gaba, " / chi.squared=", chi_da_gaba))
p_da_glut <- round(prop.test(x = c(n_da, n_glut), n = c(total_pairs, total_pairs), correct=TRUE)$p.value, 3)
chi_da_glut <- round(prop.test(x = c(n_da, n_glut), n = c(total_pairs, total_pairs), correct=TRUE)$statistic, 3)
print(paste0("DA vs glut: p=", p_da_glut, " / chi.squared=", chi_da_glut))
p_gaba_glut <- round(prop.test(x = c(n_gaba, n_glut), n = c(total_pairs, total_pairs), correct=TRUE)$p.value, 3)
chi_gaba_glut <- round(prop.test(x = c(n_gaba, n_glut), n = c(total_pairs, total_pairs), correct=TRUE)$statistic, 3)
print(paste0("GABA vs glut: p=", p_gaba_glut, " / chi.squared=", chi_gaba_glut))

print("All done!")

