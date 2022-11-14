#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
wd <- args[1]

#####################################################################
####Generating a single expression matrix of FeatureCounts method####
####################################################################

setwd(wd)

filenames = list.files(pattern = ".tsv$")
count_files = lapply(filenames, FUN = function(x) {
  read.table(x, sep = "\t", header = T)[,7]})
expr_matrix_fc <- do.call("cbind", count_files)
colnames(expr_matrix_fc) <- sub("\\..*", "", filenames)
rownames(expr_matrix_fc) <- x$Geneid
rownames(expr_matrix_fc) <- sub("\\..*", "", rownames(expr_matrix_fc))

write.csv2(expr_matrix_fc, file = paste0(wd, "raw_expr_matrix_fc.csv"))
sprintf("Raw expression matrix is saved at %s", wd)
