#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
config <- args[1]
source(as.character(config))

library(ggplot2); library(data.table)

load(expr_matrix_fc_filtered)
conds <- fread(conds_file, select = c(1:2))

#Processing matrix

expr_matrix_fc_filtered <- expr_matrix_fc_filtered + 1
expr_matrix_fc_filtered_log <- log2(expr_matrix_fc_filtered)

pcDat = prcomp(as.matrix(t(expr_matrix_fc_filtered_log)), scale = F, center = T)
var_explained <- pcDat$sdev^2/sum(pcDat$sdev^2)

pdf(file=paste0(outputdir, "/PCA_samples.labeled.pdf"), height=7, width=9)

#Trim rownames for a cleaner appearance on a plot 

rownames(pcDat$x) <- gsub("^.*P", "P", rownames(pcDat$x))

#Creating plot

theme_set(theme_bw())

pcDat$x %>% 
  as.data.frame() %>% 
  ggplot(aes(x = PC1, y = PC2, color = factor(conds$Status))) + 
  ggrepel::geom_text_repel(aes(label = rownames(pcDat$x))) +
  geom_point(size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.text = element_text(size = 12), legend.title = element_blank()) + 
  labs(x = paste0("PC1: ", round(var_explained[1]*100, 1), "%"),
       y = paste0("PC2: ", round(var_explained[2]*100, 1), "%")) +
  scale_color_manual(values = c("darkorchid4", "darkorange2")) +
  guides(color = guide_legend(override.aes = list(size=4, linetype=0))) +
  stat_ellipse(level=0.95, linetype=2)

dev.off()

