#!/usr/bin/env Rscript

#########################
####Loading libraries####
#########################

args = commandArgs(trailingOnly = TRUE)
options(stingsAsFactors = FALSE)

cat("Loading libraries...")
library(ggplot2); library(dplyr); library(stringr); library(RColorBrewer); library(haven); library(readxl); library(data.table);
library(data.table); library(ggpubr); library(gridExtra); library(cowplot)

####################
####Loading data####
####################

outputdir <- "/Users/darinaobukhova/projects/Romano_project/steroid_analysis"
cat("Loading data...")
metadata_file <- "/Users/darinaobukhova/projects/Romano_project/Romano_Overview/metadata_final.sav"
conds <- fread("/Users/darinaobukhova/projects/Romano_project/Romano_Overview/Overview_sampleID.txt",
               select = c(1:2))

#######################
####Processing data####
######################

metadata <- read_sav(metadata_file) %>% as.data.frame()
conds <- conds[order(conds$Sample_Project), ] %>% as.data.frame()
conds$SampleID <- sapply(conds$SampleID, function(x) sub(0, "", x))
metadata$ID <- sapply(metadata$ID, function(x) gsub("\\(.*", "", x)) %>% 
  sapply(function (x) sub("\\s+$", "", x))
rownames(metadata) <- metadata$ID

#Removing ID column 
metadata <- metadata[, !names(metadata) %in% c("ID")]

#Selecting columns with steroid levels in the tissue
steroid_level_tissue_data <- metadata[, grep("_T_", colnames(metadata))]

#Change colnames for better readability (removing everything after '_')
colnames(steroid_level_tissue_data) <- sapply(colnames(steroid_level_tissue_data), function(x) gsub("\\_.*", "", x))

#Rename some columns for clarity
colnames(steroid_level_tissue_data)[colnames(steroid_level_tissue_data) == "x17OHPregnenolone"] <- "17OH-Pregnenolone"
colnames(steroid_level_tissue_data)[colnames(steroid_level_tissue_data) == "x17OHProgesterone"] <- "17OH-Progesterone"
colnames(steroid_level_tissue_data)[colnames(steroid_level_tissue_data) == "AndrostEnedione"] <- "Androstenedione"

#Remove Aldosterone
steroid_level_tissue_data <- steroid_level_tissue_data[, !names(steroid_level_tissue_data) %in% c("Aldosterone", "Cortisone", "Cortisol", "Corticosterone",
                                                                                                  "x11Deoxycortisol")]

#Adding column with pregnancy status
steroid_level_tissue_data$Pregnancy_Status <- metadata$Pregnant
steroid_level_tissue_data$Pregnancy_Status <- c(rep("Non-pregnant", 20), rep("Pregnant", 20))

#Adding column with subfertility
steroid_level_tissue_data$Subfertility <- metadata$PrimSubFer
steroid_level_tissue_data$Subfertility <- ifelse(metadata$PrimSubFer == 1, "Primary subfertile", "Secondary subfertile")

#Omit columns with NA values
steroid_level_tissue_data <- steroid_level_tissue_data[,!apply(is.na(steroid_level_tissue_data), 2, any)]

#Exclude columns where all values are the same (e.g., no point of comparison)
steroid_level_tissue_data <- steroid_level_tissue_data[vapply(steroid_level_tissue_data, function(y) length(unique(y)) > 1, logical(1L))]


#Descriptive statistics
descr_stats <- steroid_level_tissue_data %>% group_by(Pregnancy) %>% summarise_all(list(mean = median, std = sd, variance = var, IQR = IQR))
write.csv2(descr_stats, file=paste0(outputdir, "/steroid_descriptive_stats.csv"))

#Stats for all steroids
steroid_tissue_levels_pvals <- apply(steroid_level_tissue_data[1:(length(steroid_level_tissue_data)-2)], 2, function(x) wilcox.test(x[1:20], x[21:40], exact = FALSE, 
                                                                      paired=FALSE, p.adjust = "bonferroni")$p.value) %>% as.data.frame()
colnames(steroid_tissue_levels_pvals) <- c("p_adjust")
write.csv2(subfert_raw_pvals, file=paste0(wd, "/steroid_tissue_levels_stats.csv"))

#####################
####Plotting data####
####################

#For pregnant and non-pregnant groups

cat("Plotting...")

#A multi-page PDF file where each page is a separate violin plot
theme_set(theme_bw())

pdf(file=paste0(outputdir, "/steroids_violinplots_updated.pdf"), height=7, width=9)

for (col in colnames(steroid_level_tissue_data[, -ncol(steroid_level_tissue_data)])) {
     ggviolin(steroid_level_tissue_data, x = "Pregnancy_Status", y = col, fill = "Pregnancy_Status",
                 palette = c("darkorchid4", "darkorange2"), add = "boxplot", add.params = list(fill = "white"), xlab=FALSE, ylab="",
                 alpha = 0.5, width = 0.8, title = col) 
                 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
                 legend.text = element_text(size = 12), legend.title = element_blank()) +
                 guides(color = guide_legend(override.aes = list(size=4, linetype=0)))
}

dev.off()

#A multipanel figure

steroid_level_plots_list <- lapply(colnames(steroid_level_tissue_data[, -ncol(steroid_level_tissue_data)]), function(var) {
  ggpar(ggviolin(steroid_level_tissue_data, x = "Pregnancy_Status", y = var, fill = "Pregnancy_Status", 
                 palette = c("darkorchid4", "darkorange2"), add = "boxplot", add.params = list(fill = "white"),
                 xlab = FALSE, ylab = FALSE, alpha = 0.5, width = 0.8, title = var) +
        geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
              plot.title = element_text(face="bold"))

)})

theme_set(theme_bw())

pdf(file=paste0(outputdir, "/Figure2_steroid_levels_tissue_updated.pdf"), height=10, width=15)

gridExtra::grid.arrange(grobs = steroid_level_plots_list, ncol = 4)

dev.off()

#Estradiol plot

ylims_estradiol <- steroid_level_tissue_data %>% 
  group_by(Pregnancy_Status) %>% 
  summarise(Q1 = quantile(Estradiol, 1/4), Q3 = quantile(Estradiol, 3/4)) %>%
  ungroup() %>%
  #get lowest Q1 and highest Q3
  summarise(lowQ1 = min(Q1), highQ3 = max(Q3))

theme_set(theme_bw())

pdf(file=paste0(outputdir, "/estradiol_scaled_y_axis.pdf"), height=7, width=9)

ggpar(ggviolin(steroid_level_tissue_data, x = "Pregnancy_Status", y = "Estradiol", fill = "Pregnancy_Status", 
               palette = c("darkorchid4", "darkorange2"), add = "boxplot", add.params = list(fill = "white"),
               xlab = FALSE, ylab = FALSE, alpha = 0.5, width = 0.8, title = "Estradiol") +
        geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
              plot.title = element_text(face="bold"))) +
        scale_y_continuous(limits = c(0, 3))
dev.off()

pdf(file=paste0(outputdir, "/estradiol_zoomed_y_axis.pdf"), height=7, width=9)

ggpar(ggviolin(steroid_level_tissue_data, x = "Pregnancy_Status", y = "Estradiol", fill = "Pregnancy_Status", 
               palette = c("darkorchid4", "darkorange2"), add = "boxplot", add.params = list(fill = "white"),
               xlab = FALSE, ylab = FALSE, alpha = 0.5, width = 0.8, title = "Estradiol") +
        geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
              plot.title = element_text(face="bold"))) +
        coord_cartesian(ylim = as.numeric(ylims_estradiol))
dev.off()

#For primary and secondary subfertile

steroid_tissue_levels_subfert_pvals <- apply(steroid_level_tissue_data[1:(length(steroid_level_tissue_data)-2)], 2, 
                                             function(x) wilcox.test(x[c(4,6,8,10,11,12,15,16,19,21,24,26,28,29,30,31,32,34,36,37)], 
                                                                     x[c(1,2,3,5,7,9,13,14,17,18,20,22,23,25,27,33,35,38,39,40)], exact = FALSE, 
                                             paired=FALSE, p.adjust = "bonferroni")$p.value) %>% as.data.frame()
colnames(steroid_tissue_levels_subfert_pvals) <- c("p_adjust")
write.csv2(subfert_raw_pvals, file=paste0(wd, "/steroid_tissue_levels_stats.csv"))

#A multipanel figure

steroid_level_plots_subfert_list <- lapply(colnames(steroid_level_tissue_data_primsubfert[1:(length(steroid_level_tissue_data_primsubfert)-2)]), function(var) {
  ggpar(ggviolin(steroid_level_tissue_data, x = "Subfertility", y = var, fill = "Subfertility", 
                 palette = c("darkseagreen4", "dodgerblue4"), add = "boxplot", add.params = list(fill = "white"),
                 xlab = FALSE, ylab = FALSE, alpha = 0.5, width = 0.8, title = var) +
          geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
                plot.title = element_text(face="bold"))
        
)})

pdf(file=paste0(outputdir, "/Figure2_steroid_levels_tissue_updated_B.pdf"), height=10, width=15)

gridExtra::grid.arrange(grobs = steroid_level_plots_subfert_list, ncol = 4)

dev.off()

#Estradiol

ylims_estradiol <- steroid_level_tissue_data %>% 
  group_by(Subfertility) %>% 
  summarise(Q1 = quantile(Estradiol, 1/4), Q3 = quantile(Estradiol, 3/4)) %>%
  ungroup() %>%
  #get lowest Q1 and highest Q3
  summarise(lowQ1 = min(Q1), highQ3 = max(Q3))

theme_set(theme_bw())

pdf(file=paste0(outputdir, "/estradiol_subfert_scaled_y_axis.pdf"), height=7, width=9)

ggpar(ggviolin(steroid_level_tissue_data, x = "Subfertility", y = "Estradiol", fill = "Subfertility", 
               palette = c("darkseagreen4", "dodgerblue4"), add = "boxplot", add.params = list(fill = "white"),
               xlab = FALSE, ylab = FALSE, alpha = 0.5, width = 0.8, title = "Estradiol") +
        geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
              plot.title = element_text(face="bold"))) +
  scale_y_continuous(limits = c(0, 3))
dev.off()

pdf(file=paste0(outputdir, "/estradiol_subfert_zoomed_y_axis.pdf"), height=7, width=9)

ggpar(ggviolin(steroid_level_tissue_data, x = "Subfertility", y = "Estradiol", fill = "Subfertility", 
               palette = c("darkseagreen4", "dodgerblue4"), add = "boxplot", add.params = list(fill = "white"),
               xlab = FALSE, ylab = FALSE, alpha = 0.5, width = 0.8, title = "Estradiol") +
        geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
              plot.title = element_text(face="bold"))) +
  coord_cartesian(ylim = as.numeric(ylims_estradiol))
dev.off()

#Primary subfertile group only

steroid_level_tissue_data_primsubfert <- steroid_level_tissue_data[steroid_level_tissue_data$Subfertility == "Primary subfertile",]

#A multipanel figure

steroid_level_plots_primsubfert_list <- lapply(colnames(steroid_level_tissue_data_primsubfert[1:(length(steroid_level_tissue_data_primsubfert)-2)]), function(var) {
  ggpar(ggviolin(steroid_level_tissue_data_primsubfert, x = "Pregnancy_Status", y = var, fill = "Pregnancy_Status", 
                 palette = c("darkorchid4", "darkorange2"), add = "boxplot", add.params = list(fill = "white"),
                 xlab = FALSE, ylab = FALSE, alpha = 0.5, width = 0.8, title = var) +
          geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
                plot.title = element_text(face="bold"))
        
)})

pdf(file=paste0(outputdir, "/Figure2_steroid_levels_tissue_updated_C.pdf"), height=10, width=15)

gridExtra::grid.arrange(grobs = steroid_level_plots_primsubfert_list, ncol = 4)

dev.off()

#Estradiol

ylims_estradiol <- steroid_level_tissue_data_primsubfert %>% 
  group_by(Pregnancy_Status) %>% 
  summarise(Q1 = quantile(Estradiol, 1/4), Q3 = quantile(Estradiol, 3/4)) %>%
  ungroup() %>%
  #get lowest Q1 and highest Q3
  summarise(lowQ1 = min(Q1), highQ3 = max(Q3))

pdf(file=paste0(outputdir, "/estradiol_primsubfert_scaled_y_axis.pdf"), height=7, width=9)

ggpar(ggviolin(steroid_level_tissue_data_primsubfert, x = "Pregnancy_Status", y = "Estradiol", fill = "Pregnancy_Status", 
               palette = c("darkorchid4", "darkorange2"), add = "boxplot", add.params = list(fill = "white"),
               xlab = FALSE, ylab = FALSE, alpha = 0.5, width = 0.8, title = "Estradiol") +
        geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
              plot.title = element_text(face="bold"))) +
  scale_y_continuous(limits = c(0, 3))
dev.off()

pdf(file=paste0(outputdir, "/estradiol_primsubfert_zoomed_y_axis.pdf"), height=7, width=9)

ggpar(ggviolin(steroid_level_tissue_data_primsubfert, x = "Pregnancy_Status", y = "Estradiol", fill = "Pregnancy_Status", 
               palette = c("darkorchid4", "darkorange2"), add = "boxplot", add.params = list(fill = "white"),
               xlab = FALSE, ylab = FALSE, alpha = 0.5, width = 0.8, title = "Estradiol") +
        geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
              plot.title = element_text(face="bold"))) +
      coord_cartesian(ylim = as.numeric(ylims_estradiol))
dev.off()

