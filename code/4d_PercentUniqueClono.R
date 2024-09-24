##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#              STEP4d: Clonality analysis II             #
#                                                        #
#    Calculates and plots the number of different clono  #
#    and the % of unique TCR sequences                   #
#                                                        #
##########################################################
#                                                        #
#                    Ambre Giguelay                      #
#                                                        #
##########################################################

rm(list = ls())

# Libraries --------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(data.table)
library(ggpubr)

# Functions --------------------------------------------------------------------

# @param seurat: seurat object
# @param dir.output: directory where all the outputs will be saved
# @param analysis.ID: ID for the output files names
# @param celltypes: vector with all the cell types included in the chosen compartment 
# @param name.annotation: name of the metadata with the cell type annotation
# @param name.sample: name of the metadata with the sample annotation
# @param name.clono: name of the metadata with the clonotype ID annotation
# @param name.cond: name of the metadata to group by (boxplot)
# @param col.samples: vector of colors for the different samples
# Saves 2 boxplots: the number of unique TCR sequence, and the % of unique TCR sequences in a sample

UniqueClono = function(seurat, dir.output, analysis.ID, celltypes, name.annotation,
                       name.sample, name.clono, name.cond, col.samples){
  cells.to.keep <- names(which(seurat[[name.annotation]][,1] %in% celltypes & seurat$ClonoStatus == "Available"))
  seurat <- subset(seurat, cells = cells.to.keep)
  samples <- unique(seurat[[name.sample]])[,1]
  meta <- seurat@meta.data
  nclono <- sapply(samples, function(x) n_distinct(meta[meta[[name.sample]] == x,][[name.clono]]))
  ncells <- table(meta[[name.sample]])
  
  
  df.plot <- data.frame(Sample = samples,
                       Ncells = as.vector(ncells[samples]),
                       Nclono = nclono[samples])
  df.plot$PercUniClono <- df.plot$Nclono/df.plot$Ncells*100
  
  condition <- table(meta[[name.sample]], meta[[name.cond]])
  condition <- setNames(colnames(condition)[apply(condition, 1, function(x) which.max(x))], rownames(condition))
  df.plot$Condition <- condition[samples]
  
  pdf(paste0(dir.output, "Plots/VDJ/BoxPlot_UniqueClono_", analysis.ID, ".pdf"), width = 6, height = 4)
  print(
    ggplot(df.plot, aes(x = Condition, y = Nclono)) +
    geom_boxplot() +
    geom_jitter(shape = 16, position = position_jitter(0.2), aes(color = Sample), size = 3) +
    scale_color_manual(values = col.samples) +
    stat_compare_means(method = "wilcox.test", paired = FALSE) +
    theme_classic() +
    theme(text = element_text(colour = "black", size = 15),
          aspect.ratio = 1/1,
          axis.line = element_line(color = "black", linewidth = 0.8),
          axis.text = element_text(colour = "black", size = 15),
          axis.ticks = element_line(color = "black", linewidth = 0.8))
  )
  dev.off()

  pdf(paste0(dir.output, "Plots/VDJ/BoxPlot_PercentageUniqueClono_", analysis.ID, ".pdf"), width = 6, height = 4)
  print(
    ggplot(df.plot, aes(x = Condition, y = PercUniClono)) +
    geom_boxplot() +
    geom_jitter(shape = 16, position = position_jitter(0.2), aes(color = Sample), size = 3) +
    scale_color_manual(values = col.samples) + 
    stat_compare_means(method = "wilcox.test", paired = FALSE) +
    theme_classic() +
    theme(text = element_text(colour = "black", size = 15),
          aspect.ratio = 1/1,
          axis.line = element_line(color = "black", linewidth = 0.8),
          axis.text = element_text(colour = "black", size = 15),
          axis.ticks = element_line(color = "black", linewidth = 0.8))
    )
  dev.off()
  
}


# Parameters --------------------------------------------------------------------
dir.output <- "..." 
samples <- ...
col.samples <- ...

# Load data --------------------------------------------------------------------
seurat.all <- readRDS(paste0(dir.output, "Objects/SeuratFinal_All.rds"))
seurat.all <- subset(seurat.all, subset = Pre_Post == "Post-ICI")

UniqueClono(seurat = seurat.all,
            celltypes = "CD4 CTL",
            name.annotation = "FinalAnnotation",
            name.sample = "Sample",
            name.clono = "clonotype_id",
            name.cond = "irAE",
            col.samples = col.samples,
            analysis.ID = "Post_CD4.CTL",
            dir.output = dir.output)

UniqueClono(seurat = seurat.all,
            celltypes = "CD8 TEM",
            name.annotation = "FinalAnnotation",
            name.sample = "Sample",
            name.clono = "clonotype_id",
            name.cond = "irAE",
            col.samples = col.samples,
            analysis.ID = "Post_CD8.TEM",
            dir.output = dir.output)

UniqueClono(seurat = seurat.all,
            celltypes = c("CD4 Naive", "CD4 TCM", "CD4 TEM", "Treg", "CD4 CTL", "CD4 Proliferating"),
            name.annotation = "FinalAnnotation",
            name.sample = "Sample",
            name.clono = "clonotype_id",
            name.cond = "irAE",
            col.samples = col.samples,
            analysis.ID = "Post_CD4",
            dir.output = dir.output)

UniqueClono(seurat = seurat.all,
            celltypes = c("CD8 Naive", "CD8 TEM", "CD8 TCM", "CD8 Proliferating"),
            name.annotation = "FinalAnnotation",
            name.sample = "Sample",
            name.clono = "clonotype_id",
            name.cond = "irAE",
            col.samples = col.samples,
            analysis.ID = "Post_CD8",
            dir.output = dir.output)