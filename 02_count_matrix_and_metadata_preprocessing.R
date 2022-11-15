#!/usr/bin/env Rscript

options(scipen=999)

args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
source(as.character(config))

#########################
####Loading libraries####
#########################

cat("Loading libraries...")

library(dplyr); library(rafalib); library(sva); library(biomaRt); library(haven); library(readxl); library(data.table); library(DESeq2); library(tidyverse)

raw_expr_matrix <- read.csv2(file = paste0(datadir, "raw_expr_matrix_fc.csv"), header = T)
steroid_genes <- read_excel(path = steroid_genes_path)
metadata <- read_sav(metadata_file,
                     col_select = c("Age_at_randomisation2", "BMI", "FSmoker_E1_C1",
                                    "FAlcohol_E1_C1", "Gravidity_E1_C2", "Parity_E1_C2",
                                    "Subfert_E1_C2", "Ferdiag_E1_C2",
                                    "CEndothick_E2_C5", "ID", "PrimSubFer", "EmbryoTransferAmount", "Timing_LH",
                                    "Pregnant")) %>% as.data.frame()

conds <- fread(conds_file, select = c(1:2))

########################################
####Processing the expression matrix####
#######################################

rownames(raw_expr_matrix) <- raw_expr_matrix$X
raw_expr_matrix <- raw_expr_matrix[, !names(raw_expr_matrix) %in% c("X")]

##################
####Filtering####
#################

#Remove rows where all values are zeros
expr_matrix_filtered <- raw_expr_matrix[rowSums(raw_expr_matrix[]) > 0, ]
sprintf("%i rows were removed", nrow(raw_expr_matrix) - nrow(expr_matrix_filtered))
sprintf("Now the expression matrix has %i gene rows", nrow(expr_matrix_filtered))

#Remove rows where values are >=5 in >75% samples

expr_matrix_filtered <- expr_matrix_filtered[rowSums(expr_matrix_filtered >= 5) >= round(ncol(expr_matrix_filtered)/100 * 75),]
sprintf("%i rows were removed", nrow(raw_expr_matrix) - nrow(expr_matrix_filtered))
sprintf("Now the expression matrix has %i gene rows", nrow(expr_matrix_filtered))

#Renaming the rows for later convenient usage with biomart
rownames(expr_matrix_filtered) <- sub("\\..*", "", rownames(expr_matrix_filtered))

##################
####Annotation####
##################

#Getting gene names
mart <- biomaRt::useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')

getinfo <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", 
             "end_position", "strand", "gene_biotype", "percentage_gene_gc_content")

gene_names_fc <- biomaRt::getBM(attributes=getinfo,
                                filters="ensembl_gene_id", 
                                values=rownames(expr_matrix_filtered), 
                                mart=mart)

#Matching the order
idx <- match(gene_names_fc$ensembl_gene_id, rownames(expr_matrix_filtered))
gene_names_fc <- gene_names_fc[order(idx),]

conds <- conds[order(conds$Sample_Project), ] %>% as.data.frame()
conds$SampleID = sapply(conds$SampleID, function(x) sub(0, "", x))
metadata$ID <- sapply(metadata$ID, function(x) gsub("\\(.*", "", x)) %>% 
  sapply(function (x) sub("\\s+$", "", x))
metadata <- metadata[metadata$ID %in% conds$SampleID, ]

rownames(metadata) <- metadata$ID

#Removing ID column 

metadata <- metadata[, !names(metadata) %in% c("ID")]

idx <- match(rownames(metadata), conds$SampleID)
metadata_ord <- metadata[order(idx),]
metadata_ord$Pregnancy_outcome <- conds$Sample_Project

#Adjusting colnames of expression matrix

colnames(expr_matrix_filtered) <- sub(0, "", colnames(expr_matrix_filtered))
colnames(expr_matrix_filtered) <- sub("\\.", "-", colnames(expr_matrix_filtered))
conds <- conds[order(match(conds$SampleID, colnames(expr_matrix_filtered))), ]
print(all(conds$SampleID == colnames(expr_matrix_filtered)))

#Excluding the samples with LH day information missing

metadata_ord <- metadata[!rownames(metadata) %in% samples_to_excl,]
expr_matrix_filtered_omitted <- expr_matrix_filtered[, !colnames(expr_matrix_filtered) %in% samples_to_excl, ]
metadata_ord <- metadata_ord[order(match(rownames(metadata_ord), colnames(expr_matrix_filtered_omitted))), ]  

#Correcting for LH day
batch <- na.omit(metadata_ord$Timing_LH)
group <- as.factor(metadata_ord$Pregnant)

corrected_data <- ComBat_seq(as.matrix(expr_matrix_filtered_omitted), batch=batch, group=group)

#Saving filtered expression matrices and corrected expression matrix

save(expr_matrix_filtered, file = paste0(datadir, "expr_matrix_filtered.RData"))
save(expr_matrix_filtered_omitted, file = paste0(datadir, "expr_matrix_filtered_samp_omitted.RData"))
save(corrected_data, gene_names_fc, file = paste0(datadir, "expr_matrix_LH_cor.RData"))
save(metadata, metadata_ord, file = paste0(datadir, "processed_metadata.RData"))

sprintf("Preprocessed data is saved at %s", datadir)
