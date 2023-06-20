print("Loading packages")
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("ggplot2")) 
suppressPackageStartupMessages(library("grid"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("dplyr"))

print("Reading in 10X filtered feature barcode matrices")
MCVTA1_data <- Read10X(data.dir="../../../raw_data/MCVTA_1/1978__filtered_feature_bc_matrix/")
MCVTA2_data <- Read10X(data.dir="../../../raw_data/MCVTA_2/1979__filtered_feature_bc_matrix/")
MCVTA3_data <- Read10X(data.dir="../../../raw_data/MCVTA_3/1986__filtered_feature_bc_matrix/")
FCVTA1_data <- Read10X(data.dir="../../../raw_data/FCVTA_1/1974__filtered_feature_bc_matrix/")
FCVTA2_data <- Read10X(data.dir="../../../raw_data/FCVTA_2/1981__filtered_feature_bc_matrix/")
FCVTA3_data <- Read10X(data.dir="../../../raw_data/FCVTA_3/1984__filtered_feature_bc_matrix/")
hungryM_data <- Read10X(data.dir="../../../raw_data/hungryM/2464__filtered_feature_bc_matrix/")
hungryF_data <- Read10X(data.dir="../../../raw_data/hungryF/2462__filtered_feature_bc_matrix/")
satedM_data <- Read10X(data.dir="../../../raw_data/satedM/2465__filtered_feature_bc_matrix/")
satedF_data <- Read10X(data.dir="../../../raw_data/satedF/2463__filtered_feature_bc_matrix/")

print("Creating Seurat objects")
MCVTA1 <- CreateSeuratObject(counts=MCVTA1_data, min.cells=10, min.features=1)
MCVTA2 <- CreateSeuratObject(counts=MCVTA2_data, min.cells=10, min.features=1)
MCVTA3 <- CreateSeuratObject(counts=MCVTA3_data, min.cells=10, min.features=1)
FCVTA1 <- CreateSeuratObject(counts=FCVTA1_data, min.cells=10, min.features=1)
FCVTA2 <- CreateSeuratObject(counts=FCVTA2_data, min.cells=10, min.features=1)
FCVTA3 <- CreateSeuratObject(counts=FCVTA3_data, min.cells=10, min.features=1)
hungryM <- CreateSeuratObject(counts=hungryM_data, min.cells=10, min.features=1)
hungryF <- CreateSeuratObject(counts=hungryF_data, min.cells=10, min.features=1)
satedM <- CreateSeuratObject(counts=satedM_data, min.cells=10, min.features=1)
satedF <- CreateSeuratObject(counts=satedF_data, min.cells=10, min.features=1)

sample_names <- c("MCVTA1", "MCVTA2", "MCVTA3", "FCVTA1", "FCVTA2", "FCVTA3", 
                  "hungryM", "hungryF", "satedM", "satedF")
samples <- c(MCVTA1, MCVTA2, MCVTA3, FCVTA1, FCVTA2, FCVTA3, hungryM, hungryF, satedM, satedF)

print("Adding metadata to Seurat objects")
for (i in 1:length(samples)) {
    
    # Add sample ID as metadata column
    samples[[i]]$sampleID <- sample_names[[i]]
    
    # Add sex as metadata column
    if (grepl("M", sample_names[[i]]) == TRUE) {
        samples[[i]]$sex <- "M"
    } else { 
    samples[[i]]$sex <- "F"
    }

    # Add hunger state as metadata column
    if (grepl("hungry", sample_names[[i]]) == TRUE) {
        samples[[i]]$condition <- "hungry"
    } else if (grepl("sated", sample_names[[i]]) == TRUE) {
        samples[[i]]$condition <- "sated"
    } else {
        samples[[i]]$condition <- "control"
    }
    
    # Add % of genes that are mitochondrial and ribosomal as metadata column
    samples[[i]][["percent.mt"]] <- PercentageFeatureSet(samples[[i]], pattern="mt-")
    samples[[i]][["percent.Rps"]] <- PercentageFeatureSet(samples[[i]], pattern="Rps")
    samples[[i]][["percent.Rpl"]] <- PercentageFeatureSet(samples[[i]], pattern="Rpl")
}

pre_filt_counts <- c()
for (i in 1:length(samples)) {
    pre <- ncol(samples[[i]][["RNA"]]@data)
    pre_filt_counts <- c(pre_filt_counts, pre)
    print(paste(sample_names[i],"has", pre, "nuclei before filtering"))
}

# Put meta data from all samples together into one df
pre_filt_meta_df <- samples[[1]]@meta.data

for (i in 2:length(samples)) {  
    df <- samples[[i]]@meta.data
    pre_filt_meta_df <- rbind(pre_filt_meta_df, df)
    }
write.csv(pre_filt_meta_df, file="../../../analysis/hungry_sated_pre_filt_metadata.csv")

print("Filtering out low quality nuclei and removing mitochondrial and ribosomal genes")
for (i in 1:length(samples)) {

    # Filter out samples that have less than 300 and more than 2500 (likely doublets) genes
    # and those that have > 5% mitochondrial genes
    samples[[i]] <- subset(x=samples[[i]], subset=nFeature_RNA>300 & nFeature_RNA<2500 & percent.mt<5)

    # Remove mitochondrial and ribosomal genes 
    all_genes <- rownames(samples[[i]])
    genes_to_remove <- rownames(samples[[i]])[startsWith(rownames(samples[[i]]), c("mt-"))]
    genes_to_remove <- c(genes_to_remove, 
                         rownames(samples[[i]])[startsWith(rownames(samples[[i]]), c("Rpl"))])
    genes_to_remove <- c(genes_to_remove, 
                         rownames(samples[[i]])[startsWith(rownames(samples[[i]]), c("Rps"))])
    genes_to_keep <- all_genes[!(all_genes %in% genes_to_remove)]
    samples[[i]] <- subset(samples[[i]], features = genes_to_keep)
}

post_filt_counts <- c()
for (i in 1:length(samples)) {
    post <- ncol(samples[[i]][["RNA"]]@data)
    post_filt_counts <- c(post_filt_counts, post)
    print(paste(sample_names[i],"has", post, "nuclei after filtering"))
    print(paste((pre_filt_counts[i]-post_filt_counts[i]), "nuclei removed"))
}

# Put meta data from all samples together into one df
post_filt_meta_df <- samples[[1]]@meta.data
for (i in 2:length(samples)) {  
    df <- samples[[i]]@meta.data
    post_filt_meta_df <- rbind(post_filt_meta_df, df)
    }
write.csv(post_filt_meta_df, file="../../../analysis/hungry_sated_post_filt_metadata.csv")

pre_counts_df <- data.frame("sample"=sample_names, 
							"nuclei"=pre_filt_counts,
							"filter"=rep("pre",length(sample_names)))
post_counts_df <- data.frame("sample"=sample_names, 
							 "nuclei"=post_filt_counts,
							 "filter"=rep("post",length(sample_names)))
counts_df <- rbind(pre_counts_df, post_counts_df)

print("Saving each individual seurat object")
for (i in 1:length(samples)) {
    savename <- paste0("../../../seurat_objects/", sample_names[i], "_preprocessed.RDS")
    saveRDS(samples[[i]], savename)
}


























