##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#      STEP4c: Clonality analysis with scRepertoire      #
#                                                        #
#    Runs scRepertoire to generate a serie of plots      #
#    (clonal diversity analysis)                         #          
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
library(scRepertoire)
library(ggpubr)

# Functions --------------------------------------------------------------------


# @param tcr.obj: TCR table (saved in 4a)
# @param dir.output: directory where all the outputs will be saved
# @param analysis.ID: ID for the output files names
# @param name.cluster: name of the chosen compartment (for outputs)
# @param col.samples: vector of colors for the different samples
# @param name.cond: name of the metadata to group by (boxplot)

RepertoireAnalysis = function(tcr.obj, dir.output, analysis.ID, name.cluster, col.samples, name.cond){
  
  # Clonotype analysis only --------------------------------
  ## quantify clonotypes 
  p1 <- quantContig(tcr.obj, cloneCall = "gene+nt", scale = T) + #Proportion of unique clonotypes (low --> clonal expansion)
    scale_fill_manual(values = col.samples) +
    theme(text = element_text(colour = "black", size = 15),
          axis.line = element_line(color = "black", linewidth = 0.8),
          axis.text = element_text(colour = "black", size = 15, angle = 30, hjust = 1),
          axis.ticks = element_line(color = "black", linewidth = 0.8))
  ggsave(filename = paste0(dir.output, "Plots/VDJ/PercentUniqueClonotype_", name.cluster, "_", analysis.ID, ".pdf"), plot = p1, width = 12, height = 4)
  
  p2 <- quantContig(tcr.obj, cloneCall = "gene+nt", scale = F) + #Absolute number of unique clontypes
    scale_fill_manual(values = col.samples) +
    theme(text = element_text(colour = "black", size = 15),
          axis.line = element_line(color = "black", linewidth = 0.8),
          axis.text = element_text(colour = "black", size = 15, angle = 30, hjust = 1),
          axis.ticks = element_line(color = "black", linewidth = 0.8))
  ggsave(filename = paste0(dir.output, "Plots/VDJ/AbundanceClonotype_", name.cluster, "_", analysis.ID, ".pdf"), plot = p2, width = 12, height = 4)
  
  ## Clonal space homeostasis
  p3 <- clonalHomeostasis(tcr.obj, cloneCall = "gene",
                         cloneTypes = c(Rare = 1e-04,
                                        Small = 0.001,
                                        Medium = 0.01,
                                        Large = 0.1,
                                        Hyperexpanded = 1)) +
    scale_fill_manual(values = c("#003366", "#006699", "#33CCCC", "#99FFFF", "#99CC99"))+
    theme(text = element_text(colour = "black", size = 15),
          axis.line = element_line(color = "black", linewidth = 0.8),
          axis.text = element_text(colour = "black", size = 15, angle = 30, hjust = 1),
          axis.ticks = element_line(color = "black", linewidth = 0.8))
  ggsave(filename = paste0(dir.output, "Plots/VDJ/ClonalHomeostasis_", name.cluster, "_", analysis.ID, ".pdf"), plot = p3, width = 12, height = 4)
  
  ## Clonal proportion
  p4 <- clonalProportion(tcr.obj, cloneCall = "gene",
                        split = c(10, 100, 1000, 10000, 30000, 1e+05)) +
    scale_fill_manual(values = c("#003366", "#006699", "#33CCCC", "#99FFFF", "#99CC99", "#CCFF99"))+
    theme(text = element_text(colour = "black", size = 15),
          axis.line = element_line(color = "black", linewidth = 0.8),
          axis.text = element_text(colour = "black", size = 15, angle = 30, hjust = 1),
          axis.ticks = element_line(color = "black", linewidth = 0.8))
  ggsave(filename = paste0(dir.output, "Plots/VDJ/ClonalProportion_", name.cluster, "_", analysis.ID, ".pdf"), plot = p4, width = 12, height = 4)
  
  ## overlap analysis
  p5 <- clonalOverlap(tcr.obj,
                     cloneCall = "gene+nt",
                     method = "morisita") +
    scale_fill_gradient(low = "lavender", high = "#663366", na.value = "white")+
    theme(text = element_text(colour = "black", size = 15),
          axis.line = element_line(color = "black", linewidth = 0.8),
          axis.text = element_text(colour = "black", size = 15, angle = 30, hjust = 1),
          axis.ticks = element_line(color = "black", linewidth = 0.8))
  ggsave(filename = paste0(dir.output, "Plots/VDJ/OverlapAnalysis_", name.cluster, "_", analysis.ID, ".pdf"), plot = p5, width = 12, height = 4)
  
  ## diversity analysis
  p6 <- clonalDiversity(tcr.obj,
                       cloneCall = "gene",
                       group.by = "sample",
                       x.axis = name.cond,
                       n.boots = 100) +
    scale_colour_manual(values = setNames(colors[1:length(tcr.obj)], names(tcr.obj))) + 
    stat_compare_means(method = "wilcox.test", paired = FALSE) 
  
  ggsave(filename = paste0(dir.output, "Plots/VDJ/DiversityAnalysis_", name.cluster, "_", analysis.ID, ".pdf"), plot = p6, width = 12, height = 4)
  
  
}

# Parameters --------------------------------------------------------------------
dir.output <- "..." 
path.cr <- "..." # path for cellranger outputs
samples <- ...
col.samples <- ...

# Analysis for Post samples only -----------------------------------------------

contig_list <- lapply(samples, function(x) fread(paste0(path.cr, x,"/outs/per_sample_outs/", x, "/vdj_t/filtered_contig_annotations.csv"),
                                                 header = T, data.table = F) )

combined <- combineTCR(contig_list, # combines the different files
                      samples = names(condition), 
                      cells ="T-AB")

for(x in names(combined)){ # Add metadata
  combined[[x]]["irAE"] <- condition[x]
}

seurat.all <- readRDS(paste0(dir.output, "Objects/SeuratFinal_All.rds"))
seurat.all <- subset(seurat.all, subset = Pre_Post == "Post-ICI")

cd8 <- rownames(filter(seurat.all@meta.data, FinalAnnotation %in% c("CD8 Naive", "CD8 TEM", "CD8 TCM", "CD8 Proliferating")))
cd4 <- rownames(filter(seurat.all@meta.data, FinalAnnotation %in% c("CD4 Naive", "CD4 TCM", "CD4 TEM", "Treg", "CD4 CTL", "CD4 Proliferating")))
cd8tem <- (filter(seurat.all@meta.data, FinalAnnotation %in% c("CD8 TEM")))
cd4ctl <- rownames(filter(seurat.all@meta.data, FinalAnnotation %in% c("CD4 CTL")))


combined.cd8 <- lapply(combined, function(x) filter(x, barcode %in% cd8))
combined.cd4 <- lapply(combined, function(x) filter(x, barcode %in% cd4))
combined.cd8tem <- lapply(combined, function(x) filter(x, barcode %in% cd8tem))
combined.cd4ctl <- lapply(combined, function(x) filter(x, barcode %in% cd4ctl))

RepertoireAnalysis(tcr.obj = combined.cd8, 
                   dir.output = dir.output, 
                   analysis.ID = "Post", 
                   name.cluster = "CD8", 
                   col.samples = col.samples,
                   name.cond = "irAE")

RepertoireAnalysis(tcr.obj = combined.cd4, 
                   dir.output = dir.output, 
                   analysis.ID = "Post", 
                   name.cluster = "CD4", 
                   col.samples = col.samples,
                   name.cond = "irAE")

RepertoireAnalysis(tcr.obj = combined.cd4ctl, 
                   dir.output = dir.output, 
                   analysis.ID = "Post", 
                   name.cluster = "CD4_CTL", 
                   col.samples = col.samples,
                   name.cond = "irAE")

RepertoireAnalysis(tcr.obj = combined.cd8tem, 
                   dir.output = dir.output, 
                   analysis.ID = "Post", 
                   name.cluster = "CD8_TEM", 
                   col.samples = col.samples,
                   name.cond = "irAE")
