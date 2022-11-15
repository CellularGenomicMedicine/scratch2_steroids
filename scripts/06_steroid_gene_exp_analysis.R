#!/usr/bin/env Rscript

options(scipen=999)
args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
source(as.character(config))

library(ggplot2); library(dplyr); library(RColorBrewer); library(PCAtools); library(VennDiagram);
library(rafalib); library(sva); library(biomaRt); library(haven); library(readxl); library(data.table); library(DESeq2); library(tidyverse)

load(corrected_data)

#########################
####Steroid presence####
########################

#Presence of steroid-metabolizing genes in the RNA-seq

steroid_genes_present <- steroid_genes[steroid_genes$Ensembl_id %in% rownames(corrected_data),]
write.csv2(steroid_genes_present, file = paste0(wd, "steroid_genes_present.csv"))
steroid_genes_absent <- steroid_genes[!steroid_genes$Ensembl_id %in% rownames(corrected_data),]
write.csv2(steroid_genes_absent, file = paste0(wd, "steroid_genes_absent.csv"))

steroid_exp <- corrected_data[rownames(corrected_data) %in% steroid_genes_present$Ensembl_id,]
steroid_exp_mat <- as.matrix(steroid_exp, ,drop = F)

idx <- match(rownames(steroid_exp), steroid_genes_present$Ensembl_id)
steroid_exp <- steroid_exp[order(idx),]
rownames(steroid_exp) <- steroid_genes_present$Gene_name
write.csv2(steroid_exp, file = paste0(wd, "/raw_steroid_exp.csv"))

transfr_steroid_exp <- as.data.frame(steroid_exp) %>% rownames_to_column() %>% gather(variable, value, -rowname) %>% spread(rowname, value)

conds_f <- conds[conds$SampleID %in% colnames(corrected_data),]
transfr_steroid_exp$variable[conds_f$Sample_Project == 'Pregnant'] <- rep("Pregnant", 17)
transfr_steroid_exp$variable[conds_f$Sample_Project == 'Non-pregnant'] <- rep("Non-Pregnant", 14)

transfr_steroid_exp <- gather(transfr_steroid_exp, AKR1C1:SULT2B1, key = "gene", value="expression")
transfr_steroid_exp$variable <- as.factor(transfr_steroid_exp$variable)
transfr_steroid_exp$gene <- as.factor(transfr_steroid_exp$gene)

pdf(file=paste0(wd, "/Figure3_preg_nonpreg_rawcounts_divided.pdf"), height=15, width=15)
ggplot(transfr_steroid_exp, aes(x = gene, y = expression, color = variable)) +
  geom_boxplot(alpha=0.8, width=0.5) +
  scale_color_manual(values = c("darkorchid4", "darkorange2")) +
  geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.title=element_blank()) +
  xlab("") +
  ylab("Expression (raw counts)") +
  theme(axis.text.x  = element_text(angle=90, vjust=0, hjust=1, size=10, colour = "black"))+
  theme(axis.text.y  = element_text(angle=0, vjust=0.5, hjust=0.5, size=10, colour = "black"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black", fill = NA, size = 1),
        panel.background = element_blank()) +
  facet_wrap(~gene, scales = "free")
dev.off()

#Wilcox test to compare steroid-metabolizing genes counts

steroid_exp_gr <- steroid_exp
colnames(steroid_exp_gr) <- ifelse(metadata_ord$Pregnant=='Yes', 'Pregnant', 'Non-pregnant')

#Order columns alphabetically
steroid_exp_gr <- steroid_exp_gr[, order(colnames(steroid_exp_gr))]
preg_nonpreg_raw_pvals <- apply(steroid_exp_gr, 1, function(x) wilcox.test(x[1:14], x[15:34], exact = FALSE, 
                                     paired=FALSE, p.adjust = "bonferroni")$p.value) %>% as.data.frame()
colnames(preg_nonpreg_raw_pvals) <- c("p_adjust")
write.csv2(preg_nonpreg_raw_pvals, file=paste0(wd, "/steroidmetgenes_preg_nonpreg_rawcounts.csv"))

#Variance Stabilizing Transformation

steroid_exp_vs <- varianceStabilizingTransformation(steroid_exp, blind = T)
write.csv2(steroid_exp_vs, file = paste0(wd, "/vst_steroid_exp.csv"))

transfr_steroid_exp_vs <- as.data.frame(steroid_exp_vs) %>% 
  rownames_to_column() %>% gather(variable, value, -rowname) %>% spread(rowname, value)

conds_f <- conds[conds$SampleID %in% colnames(corrected_data),]
transfr_steroid_exp_vs$variable[conds_f$Sample_Project == 'Pregnant'] <- rep("Pregnant", 17)
transfr_steroid_exp_vs$variable[conds_f$Sample_Project == 'Non-pregnant'] <- rep("Non-Pregnant", 14)

transfr_steroid_exp_vs <- gather(transfr_steroid_exp_vs, AKR1C1:SULT2B1, key = "gene", value="expression")
transfr_steroid_exp_vs$variable <- as.factor(transfr_steroid_exp_vs$variable)
transfr_steroid_exp_vs$gene <- as.factor(transfr_steroid_exp_vs$gene)

pdf(file=paste0(wd, "/preg_nonpreg_vst_counts.pdf"), height=10, width=12)
ggplot(transfr_steroid_exp_vs, aes(x = gene, y = expression, color = variable)) +
  geom_boxplot(alpha=0.8, width=0.5) +
  scale_color_manual(values = c("darkorchid4", "darkorange2")) +
  theme(legend.title=element_blank()) +
  xlab("") +
  ylab("Expression (vst counts)") +
  theme(axis.text.x  = element_text(angle=90, vjust=0, hjust=1, size=10, colour = "black"))+
  theme(axis.text.y  = element_text(angle=0, vjust=0.5, hjust=0.5, size=10, colour = "black"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black", fill = NA, size = 1),
        panel.background = element_blank())
dev.off()

#Steroid metabolizing genes between primary and secondary subfertile

subfert_steroid_exp <- as.data.frame(steroid_exp) %>% 
  rownames_to_column() %>% gather(variable, value, -rowname) %>% spread(rowname, value)
subfert_steroid_exp$variable <- ifelse(metadata_ord$PrimSubFer==1, 'Primary subfertile', 'Secondary subfertile')

subfert_steroid_exp <- gather(subfert_steroid_exp, AKR1C1:SULT2B1, key = "gene", value="expression")
subfert_steroid_exp$variable <- as.factor(subfert_steroid_exp$variable)
subfert_steroid_exp$gene <- as.factor(subfert_steroid_exp$gene)


pdf(file=paste0(wd, "/SuppFig1A_subfert_raw_counts_divided.pdf"), height=15, width=15)
ggplot(subfert_steroid_exp, aes(x = gene, y = expression, color = variable)) +
  geom_boxplot(alpha=0.8, width=0.3) +
  scale_color_manual(values = c("darkseagreen4", "dodgerblue4")) +
  geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.title=element_blank()) +
  xlab("") +
  ylab("Expression (raw counts)") +
  theme(axis.text.x  = element_text(angle=90, vjust=0, hjust=1, size=10, colour = "black"))+
  theme(axis.text.y  = element_text(angle=0, vjust=0.5, hjust=0.5, size=10, colour = "black"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black", fill = NA, size = 1),
        panel.background = element_blank()) +
  facet_wrap(~gene, scales = "free")
dev.off()

subfert_steroid_exp_vs <- as.data.frame(steroid_exp_vs) %>% 
  rownames_to_column() %>% gather(variable, value, -rowname) %>% spread(rowname, value)
subfert_steroid_exp_vs$variable <- ifelse(metadata_ord$PrimSubFer==1, 'Primary subfertile', 'Secondary subfertile')

subfert_steroid_exp_vs <- gather(subfert_steroid_exp_vs, AKR1C1:SULT2B1, key = "gene", value="expression")
subfert_steroid_exp_vs$variable <- as.factor(subfert_steroid_exp_vs$variable)
subfert_steroid_exp_vs$gene <- as.factor(subfert_steroid_exp_vs$gene)

pdf(file=paste0(wd, "/subfert_vst_counts.pdf"), height=8, width=12)
ggplot(subfert_steroid_exp_vs, aes(x = gene, y = expression, color = variable)) +
  geom_boxplot(alpha=0.8, width=0.3) +
  scale_color_manual(values = c("darkseagreen4", "dodgerblue4")) +
  theme(legend.title=element_blank()) +
  xlab("") +
  ylab("Expression (vst counts)") +
  theme(axis.text.x  = element_text(angle=90, vjust=0, hjust=1, size=10, colour = "black"))+
  theme(axis.text.y  = element_text(angle=0, vjust=0.5, hjust=0.5, size=10, colour = "black"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black", fill = NA, size = 1),
        panel.background = element_blank())
dev.off()
                                         
#Wilcox test to compare steroid-metabolizing genes counts

steroid_exp_gr <- steroid_exp
colnames(steroid_exp_gr) <- ifelse(metadata_ord$PrimSubFer==1, 'Primary subfertile', 'Secondary subfertile')

#Order columns alphabetically
steroid_exp_gr <- steroid_exp_gr[, order(colnames(steroid_exp_gr))]
subfert_raw_pvals <- apply(steroid_exp_gr, 1, function(x) wilcox.test(x[1:14], x[15:34], exact = FALSE, 
                                                                           paired=FALSE, p.adjust = "bonferroni")$p.value) %>% as.data.frame()
colnames(subfert_raw_pvals) <- c("p_adjust")
write.csv2(subfert_raw_pvals, file=paste0(wd, "/steroidmetgenes_subfert_prim_sec_rawcounts.csv"))

#Considering only primary subfertile group
steroid_exp_prim_subfert <- steroid_exp[, colnames(steroid_exp)[metadata_ord$PrimSubFer=="1"]]
metadata_ord_subfert <- metadata_ord[metadata_ord$PrimSubFer=="1",]
colnames(steroid_exp_prim_subfert) <- ifelse(metadata_ord_subfert$Pregnant=='Yes', 'Pregnant', 'Non-pregnant')
steroid_exp_prim_subfert <- steroid_exp_prim_subfert[, order(colnames(steroid_exp_prim_subfert))]
steroid_exp_prim_subfert_pvals <- apply(steroid_exp_prim_subfert, 1, function(x) wilcox.test(x[1:5], x[6:14], exact = FALSE, 
                                                                      paired=FALSE, p.adjust = "bonferroni")$p.value) %>% as.data.frame()
colnames(steroid_exp_prim_subfert_pvals) <- c("p_adjust")

steroid_exp_prim_subfert <- steroid_exp[, colnames(steroid_exp)[metadata_ord$PrimSubFer=="1"]] %>% as.data.frame()
steroid_exp_prim_subfert <- steroid_exp_prim_subfert %>% rownames_to_column() %>% gather(variable, value, -rowname) %>% spread(rowname, value)
steroid_exp_prim_subfert$variable <- ifelse(metadata_ord_subfert$Pregnant == 'Yes', 'Pregnant', 'Non-pregnant')
steroid_exp_prim_subfert <- gather(steroid_exp_prim_subfert, AKR1C1:SULT2B1, key = "gene", value="expression")
steroid_exp_prim_subfert$variable <- as.factor(steroid_exp_prim_subfert$variable)
steroid_exp_prim_subfert$gene <- as.factor(steroid_exp_prim_subfert$gene)

pdf(file=paste0(wd, "/SuppFig1B_primarysubfert_pregnonpregraw_counts_divided.pdf"), height=15, width=15)
ggplot(steroid_exp_prim_subfert, aes(x = gene, y = expression, color = variable)) +
  geom_boxplot(alpha=0.8, width=0.3) +
  scale_color_manual(values = c("darkorchid4", "darkorange2")) +
  geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.title=element_blank()) +
  xlab("") +
  ylab("Expression (raw counts)") +
  theme(axis.text.x  = element_text(angle=90, vjust=0, hjust=1, size=10, colour = "black"))+
  theme(axis.text.y  = element_text(angle=0, vjust=0.5, hjust=0.5, size=10, colour = "black"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black", fill = NA, size = 1),
        panel.background = element_blank()) +
  facet_wrap(~gene, scales = "free")
dev.off()

