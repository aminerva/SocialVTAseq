
calcHelper <- function(seurat_object,gene) {
    
    # INPUTS: 
        # - seurat_object: seurat object
        # - gene: ONE gene of interest; make loop if want multiple genes
    # OUTPUTS: 
        # - List containing: 
            # (1) gene name, 
            # (2) count of nuclei in the seurat object expressing the gene
            # (3) percentage of nuclei in the seurat object expressing the gene
            # (4) average expression of the gene among all nuclei in the object 
            # (5) standard error of the gene among all nuclei in the object 
            # (6) average expression of the gene in only nuclei that express the gene at any level (ie not non-expressing nuclei)
            # (7) standard error of the gene in only nuclei that express the gene at any level (ie not non-expressing nuclei)

    expr_matrix = seurat_object[['RNA']]@data
    ncells = ncol(expr_matrix)

    if(gene %in% row.names(expr_matrix)) {
        
        exp_bool <- expr_matrix[gene,]>0
        expressed_cell_names <- names(exp_bool[exp_bool==TRUE])
        count <- length(exp_bool[exp_bool==TRUE]) 
        pct <- count/ncells 
        all_nuclei_mean_expression <- mean(expr_matrix[gene,]) 
        all_nuclei_sem_expression <- std_error(expr_matrix[gene,]) 
        expressing_nuclei_mean_expression <- mean(expr_matrix[gene,expressed_cell_names]) 
        expressing_nuclei_sem_expression <- std_error(expr_matrix[gene,expressed_cell_names]) 
    
    } else {
        count <- 0 
        pct <- 0 
        all_nuclei_mean_expression <- 0 
        all_nuclei_sem_expression <- 0 
        expressing_nuclei_mean_expression <- 0 
        expressing_nuclei_sem_expression <- 0 
    }

    results <- list(gene, count, pct, round(all_nuclei_mean_expression,3), round(all_nuclei_sem_expression,3),
                    round(expressing_nuclei_mean_expression,3), round(expressing_nuclei_sem_expression,3))
    return(results)    
}



calcGeneExpressionMetrics <- function(celltype_names, genes_of_interest) {
    
    ### Function that makes a df containing important expression info from the calcHelper function for all nuclei across different celltypes ###

    # INPUTS: 
        # - celltype_names: a list of the strings of the name of the celltypes in your associated seurat objects
        # - genes_of_interest: a list of genes you want to calculate expression metrics for
    # OUTPUTS: 
        # - Dataframe containing: 
            # (1) gene names
            # (2) the number of cells expressing each gene
            # (3) the percentage of cells expressing each gene
            # (4) average expression of that gene across all cells
            # (5) standard error of expression across all cells
            # (6) average expression of that gene in only the cells that express it
            # (7) standard error of expression in only the cells that express it

    # Make an empty dataframe that you will add to
    gene_df <- data.frame("celltype"=NA, "gene"=NA, "cell_count"=NA, "pct_cells"=NA, "all_avg_exp"=NA, 
                          "all_sem_exp"=NA, "exp_avg_exp"=NA, "exp_sem_exp"=NA) 
    
    for (i in 1:length(celltype_names)) { # Iterate through all of the different cell types

        type <- celltype_names[i]
        genes <- c()
        cell_names <- c()
        cell_counts <- c()
        cell_pcts <- c()
        all_nuclei_mean_expression <- c()
        all_nuclei_sem_expression <- c()
        expressing_nuclei_mean_expression <- c()
        expressing_nuclei_sem_expression <- c()

        for (j in 1:length(genes_of_interest)) { # Iterate through all of the genes of interest

            # Calculate important metrics for a given gene in a given celltype using the calcHelper fx 
            results <- calcHelper(get(celltype_names[i]), genes_of_interest[j])
            genes <- c(genes, results[1])
            cell_counts <- c(cell_counts, results[2])
            cell_pcts <- c(cell_pcts, results[3])
            all_nuclei_mean_expression <- c(all_nuclei_mean_expression, results[4])
            all_nuclei_sem_expression <- c(all_nuclei_sem_expression, results[5])
            expressing_nuclei_mean_expression <- c(expressing_nuclei_mean_expression, results[6])
            expressing_nuclei_sem_expression <- c(expressing_nuclei_sem_expression, results[7])
        }

        # Combine into one df
        df <- data.frame("celltype" = rep(type, length(genes_of_interest)), 
                         "gene" = unlist(genes), 
                         "cell_count" = unlist(cell_counts), 
                         "pct_cells" = unlist(cell_pcts), 
                         "all_avg_exp" = unlist(all_nuclei_mean_expression), 
                         "all_sem_exp" = unlist(all_nuclei_sem_expression),
                         "exp_avg_exp" = unlist(expressing_nuclei_mean_expression), 
                         "exp_sem_exp" = unlist(expressing_nuclei_sem_expression)) 
        # Bind results from this celltype to dataframe with all celltypes
        gene_df <- rbind(gene_df, df)   
    }
    expr_df <- gene_df
    expr_df <- expr_df[-1,]
    return(expr_df)
}



calcExpressingNuclei <- function(celltype_names, genes_of_interest) {
    
    ### Function that makes a df containing expression info from the calcHelper function for only those nuclei that express each gene across different celtypes ###

    # INPUTS: 
        # - celltype_names: a list of the strings of the name of the celltypes in your associated seurat objects
        # - genes_of_interest: a list of genes you want to calculate expression metrics for
    # OUTPUTS: 
        # - Dataframe containing: 
            # (1) nucleusID
            # (2) gene name
            # (3) expression level of that gene in that nucleus
            # (4) category the gene falls into (feeding vs social)
            # (5) the celltype of that nucleus

    full_df <- data.frame("nucleusID"=NA, "gene"=NA, "expression"=NA, "celltype"=NA) 
    
    for (i in 1:length(celltype_names)) { # Iterate through the different cell types
    
        type <- celltype_names[i]

        # Make seurat object of a given celltype that contains only your genes of interest
        subsetted_obj <- subset(get(celltype_names[i]), features=genes_of_interest)

        # Get the expression data from the subsetted seurat object, transpose, and convert from matrtix to df
        subsetted_df <- data.frame(t(subsetted_obj[["RNA"]]@data))

        # Add cell IDs as a column
        subsetted_df$nucleusID <- rownames(subsetted_df)

        # Reshape for plotting
        subsetted_df <- reshape2::melt(subsetted_df, id=c("nucleusID"), 
                                       variable.name="gene", value.name="expression")

        # Only keep rows where a cell expresses a gene
        subsetted_df <- subsetted_df[subsetted_df$expression>0,]

        # Bind results from this celltype to dataframe with all celltypes
        subsetted_df$celltype <- rep(type, nrow(subsetted_df))
        full_df <- rbind(full_df, subsetted_df)
    }
    full_df <- full_df[-1,]
    return(full_df)
}


