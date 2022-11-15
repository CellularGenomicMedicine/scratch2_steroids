#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
config <- args[1]
source(as.character(config))

#Loading libraries

library(DESeq2); library(haven) 

#Loading files

load(corrected_data)

#Checking the order of the metadata columns

metadata_ord <- metadata_ord[order(match(rownames(metadata_ord), colnames(corrected_data))), ]
print(all(rownames(metadata_ord) == colnames(expr_matrix_filtered)))

#Performing DEG analysis

dds <- DESeqDataSetFromMatrix(countData = as.matrix(corrected_data),
                                 colData = metadata_ord,
                                 design = ~Pregnant)

dds$Pregnant <- relevel(dds$Pregnant, ref="No")

dds <- DESeq(dds,
             test = "Wald",
             fitType = "parametric",
             sfType = "ratio",
             modelMatrixType = "standard")
res_groups <- as.data.frame(results(dds, pAdjustMethod = "fdr"))
res_sort <- res_groups[order(res_groups$padj), ]

res_p_sig <- res_sort %>% 
  filter(padj <= 0.05) %>% 
  arrange(desc(log2FoldChange), desc(padj))

#Repeating for primary subfertile patients only

prim_subfert <- corrected_data[, colnames(corrected_data)[metadata_ord$PrimSubFer=="1"]]
metadata_ord_subfert <- metadata_ord[metadata_ord$PrimSubFer=="1",]
metadata_ord_subfert$Pregnant <- as.factor(metadata_ord_subfert$Pregnant)
print(all(rownames(metadata_ord_subfert) == colnames(prim_subfert)))

#Performing DEG analysis

dds <- DESeqDataSetFromMatrix(countData = as.matrix(prim_subfert), 
                                 colData = metadata_ord_subfert, 
                                 design = ~Pregnant)

dds$Pregnant <- relevel(dds$Pregnant, ref="No")
dds <- DESeq(dds,
             test = "Wald", 
             fitType = "parametric", 
             sfType = "ratio",
             modelMatrixType = "standard")

res_groups <- as.data.frame(results(dds, pAdjustMethod = "fdr"))
res_sort <- res_groups[order(res_groups$padj), ]

res_p_sig_primsubfert <- res_sort %>% 
  filter(padj <= 0.05) %>% 
  arrange(desc(log2FoldChange)) #significant p-adj values 

rownames(res_p_sig_primsubfert) <- gene_names$hgnc_symbol[gene_names$ensembl_gene_id %in% rownames(res_p_sig)]
sig <- res_p_sig[res_p_sig$log2FoldChange >= 0.5 | res_p_sig$log2FoldChange <= -0.5, ]

sig_upreg <- sig[sig$log2FoldChange >= 0.5, ]
sig_downreg <- sig[sig$log2FoldChange <= -0.5, ]

save(sig, sig_upreg, sig_downreg, file = paste0(wd, "PrimInfert_DEG.RData"))
write.csv2(res_p_sig_primsubfert, paste0(outputdir, "/PrimInfert_DEG_all_sign_genes.csv"))
write.csv2(sig, paste0(outputdir, "/PrimInfert_DEG_sign_logFCthres_genes.csv"))

